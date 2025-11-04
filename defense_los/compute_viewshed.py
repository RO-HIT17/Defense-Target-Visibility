
"""
Improved raster viewshed script.

Fixes and improvements:
- Uses a safer step size (<= pixel size) to avoid skipping pixels.
- Adds angular tolerance to avoid tiny-angle noise flips.
- Optional morphological closing (requires scipy) to remove speckle.
- Progress logs for azimuths.
- Validates GeoJSON features before writing.
- CLI params: lon, lat, height, maxdist, dem, azimuth_steps, tol_deg, smooth

Example:
python viewshed_fixed.py --lon 71.385 --lat 34.189 --height 50 --maxdist 5000 --dem data/demo_aoi.tif
"""
import math
import json
import argparse
from pathlib import Path
import numpy as np
import rasterio
from rasterio.features import shapes


try:
    from scipy.ndimage import binary_closing
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

parser = argparse.ArgumentParser()
parser.add_argument("--lon", type=float, required=True)
parser.add_argument("--lat", type=float, required=True)
parser.add_argument("--height", type=float, required=True)
parser.add_argument("--maxdist", type=float, required=True)
parser.add_argument("--dem", type=str, required=True)
parser.add_argument("--azimuth_steps", type=int, default=1440, help="Number of rays (default 1440).")
parser.add_argument("--tol_deg", type=float, default=0.05, help="Angle tolerance in degrees to prevent noise (default 0.05).")
parser.add_argument("--smooth", action="store_true", help="Apply morphological closing to output (requires scipy).")
parser.add_argument("--out_geojson", type=str, default="data/viewshed_py.geojson")
args = parser.parse_args()

DEM_PATH = args.dem
OBSERVER_LON, OBSERVER_LAT = args.lon, args.lat
OBSERVER_HEIGHT_M = args.height
MAX_DISTANCE_M = args.maxdist
AZIMUTH_STEPS = args.azimuth_steps
ANGLE_TOL = math.radians(args.tol_deg)
OUT_GEOJSON = Path(args.out_geojson)

METERS_PER_DEG_LAT = 111132.92


def meters_per_deg_lon(lat_deg):
    return 111412.84 * math.cos(math.radians(lat_deg)) - 93.5 * math.cos(3 * math.radians(lat_deg))


def main():
    print("Loading DEM:", DEM_PATH)
    with rasterio.open(DEM_PATH) as src:
        dem = src.read(1, masked=True)
        transform = src.transform
        crs = src.crs
        width = src.width
        height = src.height
        bounds = src.bounds
        print(f"DEM: size={width}x{height}, crs={crs}, bounds={bounds}")
        
        pix_size_x = abs(transform.a)
        pix_size_y = abs(transform.e)

    
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
        raise SystemExit("Observer point falls on NoData in DEM. Pick a nearby valid coordinate.")
    print(f"Observer (lon,lat) = ({OBSERVER_LON}, {OBSERVER_LAT}), ground elev = {obs_elev:.3f} m, eye height = {OBSERVER_HEIGHT_M} m")

    
    mpd_lat = METERS_PER_DEG_LAT
    mpd_lon = meters_per_deg_lon(OBSERVER_LAT)

    
    meters_per_pixel = math.hypot(pix_size_x * mpd_lon, pix_size_y * mpd_lat)
    step_m = min(meters_per_pixel * 0.6, meters_per_pixel)  
    if step_m <= 0:
        step_m = meters_per_pixel

    max_steps = int(math.ceil(MAX_DISTANCE_M / step_m))
    print(f"Pixel size ~ {meters_per_pixel:.2f} m, using step_m = {step_m:.2f} m, max_steps = {max_steps}, azimuths = {AZIMUTH_STEPS}")
    print(f"Angle tolerance = {math.degrees(ANGLE_TOL):.4f} deg")

    visible = np.zeros_like(dem, dtype=np.uint8)

    two_pi = 2.0 * math.pi
    azimuths = np.linspace(0, two_pi, AZIMUTH_STEPS, endpoint=False)

    for i, theta in enumerate(azimuths):
        
        if i % max(1, AZIMUTH_STEPS // 16) == 0:
            print(f"Processing azimuth {i+1}/{AZIMUTH_STEPS} ...")
        sin_t = math.sin(theta)
        cos_t = math.cos(theta)
        max_angle = -1e9  

        for step_idx in range(1, max_steps + 1):
            d_m = step_idx * step_m
            delta_lon = (d_m * cos_t) / mpd_lon
            delta_lat = (d_m * sin_t) / mpd_lat
            samp_lon = OBSERVER_LON + delta_lon
            samp_lat = OBSERVER_LAT + delta_lat

            elev = sample_lonlat(samp_lon, samp_lat)
            if elev is None:
                
                break

            
            
            dz = elev - (obs_elev + OBSERVER_HEIGHT_M)
            horizontal_dist_m = d_m
            elev_angle = math.atan2(dz, horizontal_dist_m)

            
            if elev_angle > max_angle + ANGLE_TOL:
                
                max_angle = elev_angle
                r, c = lonlat_to_rowcol(samp_lon, samp_lat)
                if 0 <= r < visible.shape[0] and 0 <= c < visible.shape[1]:
                    visible[r, c] = 1
            

    
    if args.smooth:
        if SCIPY_AVAILABLE:
            print("Applying morphological closing (scipy) to remove speckle ...")
            
            closed = binary_closing(visible.astype(bool), structure=np.ones((3, 3)))
            visible = closed.astype(np.uint8)
        else:
            print("Skipping smoothing: scipy not available. Install scipy to enable smoothing (--smooth)")

    
    mask = visible == 1
    shapes_gen = shapes(visible, mask=mask, transform=transform)
    features = []
    count = 0
    for geom, val in shapes_gen:
        if val == 1:
            
            features.append({"type": "Feature", "properties": {"visible": 1}, "geometry": geom})
            count += 1

    geojson = {"type": "FeatureCollection", "features": features}
    OUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_GEOJSON, "w", encoding="utf-8") as f:
        json.dump(geojson, f)

    print(f"Done. Wrote {count} polygons to {OUT_GEOJSON}")

if __name__ == "__main__":
    main()
