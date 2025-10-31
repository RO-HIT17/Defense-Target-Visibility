import math
import numpy as np
import rasterio
from rasterio.features import shapes
import json
from shapely.geometry import shape, mapping
from rasterio.transform import from_origin


def meters_per_deg_lon(lat_deg):
    return 111412.84 * math.cos(math.radians(lat_deg)) - 93.5 * math.cos(3 * math.radians(lat_deg))


def run(
    DEM_PATH,
    OUT_TIF,
    OUT_GEOJSON,
    OBSERVER_LON,
    OBSERVER_LAT,
    OBSERVER_HEIGHT_M=50.0,
    MAX_DISTANCE_M=5000,
    AZIMUTH_STEPS=720,
    STEP_M=None
):
    print(f"\n🗺️ Starting viewshed calculation for:")
    print(f"  DEM: {DEM_PATH}")
    print(f"  Observer: lon={OBSERVER_LON}, lat={OBSERVER_LAT}")
    print(f"  Observer height: {OBSERVER_HEIGHT_M} m")
    print(f"  Max distance: {MAX_DISTANCE_M} m")
    print(f"  Azimuth steps: {AZIMUTH_STEPS}")
    print("-" * 80)

    with rasterio.open(DEM_PATH) as src:
        dem = src.read(1, masked=True)
        transform = src.transform
        crs = src.crs
        width = src.width
        height = src.height
        bounds = src.bounds

    print(f"✅ DEM loaded successfully.")
    print(f"  Dimensions: {width} x {height}")
    print(f"  CRS: {crs}")
    print(f"  Bounds: {bounds}")
    print(f"  Pixel size (x, y): {transform.a:.6f}, {transform.e:.6f}")
    print("-" * 80)

    inv = ~transform

    def lonlat_to_rowcol(x, y):
        col, row = ~transform * (x, y)
        return int(round(row)), int(round(col))

    def sample_lonlat(x, y):
        r, c = lonlat_to_rowcol(x, y)
        if r < 0 or r >= dem.shape[0] or c < 0 or c >= dem.shape[1]:
            return None
        val = dem[r, c]
        if np.ma.is_masked(val):
            return None
        return float(val)

    obs_elev = sample_lonlat(OBSERVER_LON, OBSERVER_LAT)
    if obs_elev is None:
        raise SystemExit("🚫 Observer point falls on NoData in DEM. Pick a nearby valid coordinate.")

    print(f"👁️ Observer ground elevation: {obs_elev:.2f} m")
    print(f"👁️ Total observer height (ground + tower): {obs_elev + OBSERVER_HEIGHT_M:.2f} m")
    print("-" * 80)

    visible = np.zeros_like(dem, dtype=np.uint8)

    # Distance calculations
    METERS_PER_DEG_LAT = 111132.92
    mpd_lat = METERS_PER_DEG_LAT
    mpd_lon = meters_per_deg_lon(OBSERVER_LAT)

    pix_size_x = transform.a
    pix_size_y = -transform.e
    meters_per_pixel = math.hypot(pix_size_x * mpd_lon, pix_size_y * mpd_lat)
    step_m = STEP_M or (meters_per_pixel * 0.9)

    print(f"📏 Meters per degree (lat): {mpd_lat:.2f}")
    print(f"📏 Meters per degree (lon): {mpd_lon:.2f}")
    print(f"📏 Estimated meters per pixel: {meters_per_pixel:.2f}")
    print(f"📏 Step distance: {step_m:.2f} m per ray step")
    print("-" * 80)

    max_steps = int(math.ceil(MAX_DISTANCE_M / step_m))
    print(f"🔁 Total steps per ray: {max_steps}")

    two_pi = 2 * math.pi
    for i, theta in enumerate(np.linspace(0, two_pi, AZIMUTH_STEPS, endpoint=False)):
        if i % (AZIMUTH_STEPS // 8) == 0:
            print(f"  Processing azimuth {i}/{AZIMUTH_STEPS} ({(i/AZIMUTH_STEPS)*100:.1f}%) ...")

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

            dz = elev - obs_elev
            elev_angle = math.atan2(dz - OBSERVER_HEIGHT_M, d_m)
            if elev_angle > max_angle + 1e-12:
                max_angle = elev_angle
                r, c = lonlat_to_rowcol(samp_lon, samp_lat)
                if 0 <= r < visible.shape[0] and 0 <= c < visible.shape[1]:
                    visible[r, c] = 1

            if max_angle > math.radians(80):
                break

    print("✅ Visibility analysis complete.")
    print("-" * 80)

    # Save raster output
    out_meta = {
        'driver': 'GTiff',
        'height': visible.shape[0],
        'width': visible.shape[1],
        'count': 1,
        'dtype': 'uint8',
        'crs': crs,
        'transform': transform,
        'compress': 'LZW'
    }
    with rasterio.open(OUT_TIF, 'w', **out_meta) as dst:
        dst.write(visible, 1)

    print(f"💾 Saved visibility raster → {OUT_TIF}")

    # Convert visible area to GeoJSON
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
    with open(OUT_GEOJSON, 'w', encoding='utf-8') as f:
        json.dump(geojson, f)

    print(f"💾 Saved GeoJSON output → {OUT_GEOJSON}")
    print("🎯 Done! Viewshed generation completed successfully.")


if __name__ == "__main__":
    print("\n=== Dynamic Viewshed Runner ===")
    DEM_PATH = input("Enter DEM path (e.g. data/demo_crop.tif): ").strip() or r"data\demo_crop.tif"
    OUT_TIF = input("Enter output TIFF path: ").strip() or r"data\viewshed_py.tif"
    OUT_GEOJSON = input("Enter output GeoJSON path: ").strip() or r"data\viewshed_py.geojson"
    OBSERVER_LON = float(input("Enter observer longitude: ").strip() or 71.385)
    OBSERVER_LAT = float(input("Enter observer latitude: ").strip() or 34.189)
    OBSERVER_HEIGHT_M = float(input("Enter observer height (m): ").strip() or 50.0)
    MAX_DISTANCE_M = float(input("Enter max distance (m): ").strip() or 5000)
    AZIMUTH_STEPS = int(input("Enter azimuth steps: ").strip() or 720)
    STEP_M = input("Enter step distance (m) or leave blank for auto: ").strip()
    STEP_M = float(STEP_M) if STEP_M else None

    run(
        DEM_PATH,
        OUT_TIF,
        OUT_GEOJSON,
        OBSERVER_LON,
        OBSERVER_LAT,
        OBSERVER_HEIGHT_M,
        MAX_DISTANCE_M,
        AZIMUTH_STEPS,
        STEP_M
    )
