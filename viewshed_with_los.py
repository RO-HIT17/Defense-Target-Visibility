import math
import numpy as np
import rasterio
from rasterio.features import shapes
import json
from shapely.geometry import shape, mapping, box
import logging
from rasterio.transform import from_origin

# --------------------------------------------
# CONFIGURATION
# --------------------------------------------
DEM_PATH = r"data\demo_aoi.tif"
OUT_TIF = r"data\viewshed_py.tif"
OUT_GEOJSON = r"data\viewshed_py.geojson"
OUT_LOS_GEOJSON = r"data\viewshed_los.geojson"

# Observer (you can change these or wire to CLI/Streamlit)
OBSERVER_LON = 71.985
OBSERVER_LAT = 34.189
OBSERVER_HEIGHT_M = 40.0

# Target for LOS check (change coordinates and optional height)
TARGET_LON = 71.983
TARGET_LAT = 34.182
TARGET_HEIGHT_M = 10.0   # height above ground for target (e.g., top of tower)

MAX_DISTANCE_M = 500
AZIMUTH_STEPS = 720
STEP_M = None

METERS_PER_DEG_LAT = 111132.92

# --------------------------------------------
# LOGGING SETUP
# --------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)

# --------------------------------------------
# HELPER FUNCTIONS
# --------------------------------------------
def meters_per_deg_lon(lat_deg):
    """Approx meters per degree longitude at given latitude (ellipsoidal approx)."""
    return 111412.84 * math.cos(math.radians(lat_deg)) - 93.5 * math.cos(3 * math.radians(lat_deg))


def run():
    logging.info("=== VIEWSHED + LOS CALCULATION STARTED ===")
    logging.info(f"DEM Path: {DEM_PATH}")
    logging.info(f"Observer: Lon={OBSERVER_LON}, Lat={OBSERVER_LAT}, Height={OBSERVER_HEIGHT_M}m")
    logging.info(f"Target:   Lon={TARGET_LON}, Lat={TARGET_LAT}, Height={TARGET_HEIGHT_M}m")

    # Open DEM (assumes DEM in lon/lat EPSG:4326; see notes if not)
    with rasterio.open(DEM_PATH) as src:
        dem = src.read(1, masked=True)
        transform = src.transform
        crs = src.crs
        bounds = src.bounds

        logging.info(f"DEM Metadata:")
        logging.info(f"  CRS: {crs}")
        logging.info(f"  Width x Height: {src.width} x {src.height}")
        logging.info(f"  Bounds: {bounds}")
        logging.info(f"  Transform: {transform}")

    # Inverse affine for lon/lat -> col/row mapping
    inv = ~transform

    def lonlat_to_rowcol(x, y):
        """Return (row, col) integer indices for lon/lat coordinate."""
        col, row = ~transform * (x, y)
        return int(round(row)), int(round(col))

    def sample_lonlat(x, y):
        """Sample DEM at lon/lat. Returns float elevation or None if outside/NoData."""
        r, c = lonlat_to_rowcol(x, y)
        if r < 0 or r >= dem.shape[0] or c < 0 or c >= dem.shape[1]:
            return None
        val = dem[r, c]
        if np.ma.is_masked(val):
            return None
        return float(val)

    # Validate observer
    obs_elev = sample_lonlat(OBSERVER_LON, OBSERVER_LAT)
    if obs_elev is None:
        logging.error("Observer point falls on NoData in DEM.")
        raise SystemExit("Invalid observer coordinate.")

    logging.info(f"Observer ground elevation: {obs_elev:.2f} m")

    # ---------- Viewshed (kept from your code) ----------
    visible = np.zeros_like(dem, dtype=np.uint8)
    mpd_lat = METERS_PER_DEG_LAT
    mpd_lon = meters_per_deg_lon(OBSERVER_LAT)
    pix_size_x = transform.a
    pix_size_y = -transform.e
    meters_per_pixel = math.hypot(pix_size_x * mpd_lon, pix_size_y * mpd_lat)
    step_m = meters_per_pixel * 0.9 if STEP_M is None else STEP_M
    max_steps = int(math.ceil(MAX_DISTANCE_M / step_m))
    logging.info(f"Computed step size: {step_m:.2f} m, Max steps: {max_steps}")

    two_pi = 2 * math.pi
    logging.info(f"Processing {AZIMUTH_STEPS} azimuth directions...")

    for i, theta in enumerate(np.linspace(0, two_pi, AZIMUTH_STEPS, endpoint=False)):
        if i % 100 == 0:
            logging.info(f"  → Azimuth {i}/{AZIMUTH_STEPS}")
        sin_t = math.sin(theta)
        cos_t = math.cos(theta)
        max_angle = -999.0

        for step_idx in range(1, max_steps + 1):
            d_m = step_idx * step_m
            delta_lon = (d_m * cos_t) / mpd_lon
            delta_lat = (d_m * sin_t) / mpd_lat
            samp_lon = OBSERVER_LON + delta_lon
            samp_lat = OBSERVER_LAT + delta_lat
            elev = sample_lonlat(samp_lon, samp_lat)
            if elev is None:
                break

            horizontal_dist_m = d_m
            dz = elev - obs_elev
            elev_angle = math.atan2(dz - OBSERVER_HEIGHT_M, horizontal_dist_m)
            if elev_angle > max_angle + 1e-12:
                max_angle = elev_angle
                r, c = lonlat_to_rowcol(samp_lon, samp_lat)
                if 0 <= r < visible.shape[0] and 0 <= c < visible.shape[1]:
                    visible[r, c] = 1
            if max_angle > math.radians(80):
                break

    # --------------------------------------------
    # WRITE OUTPUT RASTER
    # --------------------------------------------
    out_meta = {
        "driver": "GTiff",
        "height": visible.shape[0],
        "width": visible.shape[1],
        "count": 1,
        "dtype": "uint8",
        "crs": crs,
        "transform": transform,
        "compress": "LZW"
    }

    with rasterio.open(OUT_TIF, "w", **out_meta) as dst:
        dst.write(visible, 1)

    logging.info(f"Saved visibility raster to {OUT_TIF}")

    # --------------------------------------------
    # EXTRACT POLYGONS TO GEOJSON
    # --------------------------------------------
    mask = visible == 1
    shapes_gen = shapes(visible, mask=mask, transform=transform)
    features = []
    for geom, val in shapes_gen:
        if val == 1:
            features.append({
                "type": "Feature",
                "properties": {"visible": int(val)},
                "geometry": geom
            })

    geojson = {"type": "FeatureCollection", "features": features}
    with open(OUT_GEOJSON, "w", encoding="utf-8") as f:
        json.dump(geojson, f)

    logging.info(f"Saved GeoJSON with {len(features)} polygons to {OUT_GEOJSON}")

    # --------------------------------------------
    # LINE-OF-SIGHT (LOS) IMPLEMENTATION
    # --------------------------------------------
    logging.info("Running Line-Of-Sight (LOS) check...")

    def haversine_dist_m(lon1, lat1, lon2, lat2):
        # Great-circle distance approximation (meters)
        R = 6371000.0
        lon1r, lat1r, lon2r, lat2r = map(math.radians, (lon1, lat1, lon2, lat2))
        dlon = lon2r - lon1r
        dlat = lat2r - lat1r
        a = math.sin(dlat/2)**2 + math.cos(lat1r)*math.cos(lat2r)*math.sin(dlon/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        return R * c

    def sample_line(observer_lon, observer_lat, target_lon, target_lat, step_m_local=None):
        """Generate sample lon/lat coordinates along straight line from observer to target."""
        total_dist = haversine_dist_m(observer_lon, observer_lat, target_lon, target_lat)
        if total_dist == 0:
            return []
        if step_m_local is None:
            step = max(1.0, step_m)  # fallback
        else:
            step = step_m_local
        n_steps = int(math.ceil(total_dist / step))
        samples = []
        for i in range(0, n_steps + 1):
            frac = i / n_steps
            # simple linear interpolation in lon/lat — OK for short distances
            lon = observer_lon + frac * (target_lon - observer_lon)
            lat = observer_lat + frac * (target_lat - observer_lat)
            d = frac * total_dist
            samples.append((frac, d, lon, lat))
        return samples

    def los_check(observer_lon, observer_lat, observer_height_m,
                  target_lon, target_lat, target_height_m, step_m_local=None):
        """
        Returns:
          visible (bool),
          blocking_info: None if visible else dict with {lon, lat, elev, frac, distance_m}
        """
        obs_ground_z = sample_lonlat(observer_lon, observer_lat)
        tgt_ground_z = sample_lonlat(target_lon, target_lat)
        if obs_ground_z is None or tgt_ground_z is None:
            logging.warning("Observer or target falls on NoData in DEM.")
            return False, {"reason": "nodata"}

        obs_z = obs_ground_z + observer_height_m
        tgt_z = tgt_ground_z + target_height_m
        samples = sample_line(observer_lon, observer_lat, target_lon, target_lat, step_m_local=step_m_local)

        for frac, d, lon, lat in samples:
            # skip observer point (i=0) and target point evaluation (we may treat target as visible even if exactly on top)
            if d == 0:
                continue
            elev_sample = sample_lonlat(lon, lat)
            if elev_sample is None:
                # If we cross NoData, conservatively treat as blocked (or you may choose to skip)
                return False, {"reason": "nodata_cross", "lon": lon, "lat": lat}
            # elevation of the line at this fraction
            z_line = obs_z + frac * (tgt_z - obs_z)
            # if terrain at sample point is above the LOS line -> blocked
            if elev_sample > z_line + 1e-6:
                return False, {
                    "lon": lon,
                    "lat": lat,
                    "elev": elev_sample,
                    "frac": frac,
                    "distance_m": d
                }
        return True, None

    visible_los, blocking_info = los_check(
        OBSERVER_LON, OBSERVER_LAT, OBSERVER_HEIGHT_M,
        TARGET_LON, TARGET_LAT, TARGET_HEIGHT_M,
        step_m_local=step_m
    )

    logging.info(f"LOS result: visible={visible_los}, blocking_info={blocking_info}")

    # Prepare LOS GeoJSON: LineString + markers
    los_features = []
    line_feat = {
        "type": "Feature",
        "properties": {
            "type": "line_of_sight",
            "visible": bool(visible_los),
            "observer": {"lon": OBSERVER_LON, "lat": OBSERVER_LAT, "height_m": OBSERVER_HEIGHT_M},
            "target": {"lon": TARGET_LON, "lat": TARGET_LAT, "height_m": TARGET_HEIGHT_M}
        },
        "geometry": {
            "type": "LineString",
            "coordinates": [
                [OBSERVER_LON, OBSERVER_LAT],
                [TARGET_LON, TARGET_LAT]
            ]
        }
    }
    los_features.append(line_feat)

    # Observer and target markers
    los_features.append({
        "type": "Feature",
        "properties": {"type": "observer", "height_m": OBSERVER_HEIGHT_M},
        "geometry": {"type": "Point", "coordinates": [OBSERVER_LON, OBSERVER_LAT]}
    })
    los_features.append({
        "type": "Feature",
        "properties": {"type": "target", "height_m": TARGET_HEIGHT_M, "visible": bool(visible_los)},
        "geometry": {"type": "Point", "coordinates": [TARGET_LON, TARGET_LAT]}
    })

    # If blocked, add blocking point
    if not visible_los and blocking_info and "lon" in blocking_info:
        los_features.append({
            "type": "Feature",
            "properties": {"type": "blocker", "elev": blocking_info.get("elev"), "distance_m": blocking_info.get("distance_m")},
            "geometry": {"type": "Point", "coordinates": [blocking_info["lon"], blocking_info["lat"]]}
        })

    los_geojson = {"type": "FeatureCollection", "features": los_features}
    with open(OUT_LOS_GEOJSON, "w", encoding="utf-8") as f:
        json.dump(los_geojson, f)
    logging.info(f"Saved LOS GeoJSON to {OUT_LOS_GEOJSON}")

    # --------------------------------------------
    # GEOJSON BOUNDS (viewshed)
    # --------------------------------------------
    try:
        minx = miny = float("inf")
        maxx = maxy = float("-inf")
        for feat in features:
            geom = shape(feat["geometry"])
            bx, by, Bx, By = geom.bounds
            minx, miny = min(minx, bx), min(miny, by)
            maxx, maxy = max(maxx, Bx), max(maxy, By)
        logging.info(f"GeoJSON Bounds: ({minx:.6f}, {miny:.6f}) - ({maxx:.6f}, {maxy:.6f})")
    except Exception as e:
        logging.warning(f"Could not compute GeoJSON bounds: {e}")

    logging.info("=== VIEWSHED + LOS CALCULATION COMPLETED SUCCESSFULLY ===")


if __name__ == "__main__":
    run()
