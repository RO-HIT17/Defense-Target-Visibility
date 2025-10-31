#!/usr/bin/env python3
import math
import json
from pathlib import Path
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterio.warp import transform_geom

DEM_PATH = Path("data/demo_crop.tif")
OUTPUT_DIR = Path("data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OBS_POINTS = [
    {"name": "obsA", "lon": 71.485, "lat": 34.289, "height_m": 50.0},
    {"name": "obsB", "lon": 71.385, "lat": 34.189, "height_m": 10.0},
]

MAX_DISTANCE_M = 5000.0
AZIMUTH_STEPS = 720
STEP_M = None


def meters_per_deg_lon(lat_deg):
    return 111412.84 * math.cos(math.radians(lat_deg)) - 93.5 * math.cos(3 * math.radians(lat_deg))


def run_for_observer(obs, dem_path, out_tif, out_geojson, max_distance_m=5000, azimuth_steps=720, step_m=None):
    with rasterio.open(dem_path) as src:
        dem = src.read(1, masked=True)
        transform = src.transform
        crs = src.crs

    observer_x, observer_y = obs["lon"], obs["lat"]

    def xy_to_rowcol(x, y):
        col, row = ~transform * (x, y)
        return int(round(row)), int(round(col))

    def sample_xy(x, y):
        r, c = xy_to_rowcol(x, y)
        if r < 0 or r >= dem.shape[0] or c < 0 or c >= dem.shape[1]:
            return None
        val = dem[r, c]
        if np.ma.is_masked(val):
            return None
        return float(val)

    obs_elev = sample_xy(observer_x, observer_y)
    if obs_elev is None:
        raise SystemExit("Observer point falls on NoData in DEM.")

    if step_m is None:
        pix_size_x, pix_size_y = transform.a, -transform.e
        mpd_lat = 111132.92
        mpd_lon = meters_per_deg_lon(obs["lat"])
        meters_per_pixel = math.hypot(abs(pix_size_x) * mpd_lon, abs(pix_size_y) * mpd_lat)
        step_m = meters_per_pixel * 0.9

    max_steps = int(math.ceil(max_distance_m / step_m))
    visible = np.zeros_like(dem, dtype=np.uint8)
    two_pi = 2 * math.pi

    for theta in np.linspace(0, two_pi, azimuth_steps, endpoint=False):
        sin_t, cos_t = math.sin(theta), math.cos(theta)
        max_angle = -999.0
        for step_idx in range(1, max_steps + 1):
            d_m = step_idx * step_m
            delta_lon = (d_m * cos_t) / meters_per_deg_lon(obs["lat"])
            delta_lat = (d_m * sin_t) / 111132.92
            samp_x = observer_x + delta_lon
            samp_y = observer_y + delta_lat

            elev = sample_xy(samp_x, samp_y)
            if elev is None:
                break

            elev_angle = math.atan2(elev - obs_elev - obs["height_m"], d_m)
            if elev_angle > max_angle + 1e-12:
                max_angle = elev_angle
                r, c = xy_to_rowcol(samp_x, samp_y)
                if 0 <= r < visible.shape[0] and 0 <= c < visible.shape[1]:
                    visible[r, c] = 1
            if max_angle > math.radians(80):
                break

    out_meta = {
        "driver": "GTiff",
        "height": visible.shape[0],
        "width": visible.shape[1],
        "count": 1,
        "dtype": "uint8",
        "crs": crs,
        "transform": transform,
        "compress": "LZW",
    }
    with rasterio.open(out_tif, "w", **out_meta) as dst:
        dst.write(visible, 1)

    mask = visible == 1
    features = []
    for geom, val in shapes(visible, mask=mask, transform=transform):
        if val == 1:
            geom_wgs84 = transform_geom(crs, "EPSG:4326", geom, precision=6)
            features.append({"type": "Feature", "properties": {"visible": int(val)}, "geometry": geom_wgs84})

    with open(out_geojson, "w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)


def main():
    for obs in OBS_POINTS:
        name = obs["name"]
        out_tif = OUTPUT_DIR / f"viewshed_{name}.tif"
        out_geojson = OUTPUT_DIR / f"viewshed_{name}.geojson"
        run_for_observer(obs, DEM_PATH, out_tif, out_geojson, MAX_DISTANCE_M, AZIMUTH_STEPS, STEP_M)


if __name__ == "__main__":
    main()
