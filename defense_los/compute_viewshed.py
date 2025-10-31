import math, json, rasterio, argparse
import numpy as np
from rasterio.features import shapes

parser = argparse.ArgumentParser()
parser.add_argument("--lon", type=float)
parser.add_argument("--lat", type=float)
parser.add_argument("--height", type=float)
parser.add_argument("--maxdist", type=float)
parser.add_argument("--dem", type=str)
args = parser.parse_args()

DEM_PATH = args.dem
OBSERVER_LON, OBSERVER_LAT = args.lon, args.lat
OBSERVER_HEIGHT_M = args.height
MAX_DISTANCE_M = args.maxdist

METERS_PER_DEG_LAT = 111132.92

def meters_per_deg_lon(lat_deg):
    return 111412.84 * math.cos(math.radians(lat_deg)) - 93.5 * math.cos(3 * math.radians(lat_deg))

with rasterio.open(DEM_PATH) as src:
    dem = src.read(1, masked=True)
    transform = src.transform
    crs = src.crs

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
visible = np.zeros_like(dem, dtype=np.uint8)

mpd_lat = METERS_PER_DEG_LAT
mpd_lon = meters_per_deg_lon(OBSERVER_LAT)
pix_size_x, pix_size_y = transform.a, -transform.e
step_m = math.hypot(pix_size_x * mpd_lon, pix_size_y * mpd_lat)
max_steps = int(math.ceil(MAX_DISTANCE_M / step_m))

for theta in np.linspace(0, 2*math.pi, 720, endpoint=False):
    sin_t, cos_t = math.sin(theta), math.cos(theta)
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

mask = visible == 1
shapes_gen = shapes(visible, mask=mask, transform=transform)
features = [
    {"type": "Feature", "properties": {"visible": 1}, "geometry": geom}
    for geom, val in shapes_gen if val == 1
]

with open("data/viewshed_py.geojson", "w") as f:
    json.dump({"type": "FeatureCollection", "features": features}, f)
