import argparse, rasterio, math, json, numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--lon1", type=float)
parser.add_argument("--lat1", type=float)
parser.add_argument("--lon2", type=float)
parser.add_argument("--lat2", type=float)
parser.add_argument("--h1", type=float)
parser.add_argument("--h2", type=float)
parser.add_argument("--dem", type=str)
args = parser.parse_args()

with rasterio.open(args.dem) as src:
    dem = src.read(1, masked=True)
    transform = src.transform

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

N = 200
lons = np.linspace(args.lon1, args.lon2, N)
lats = np.linspace(args.lat1, args.lat2, N)
dists = np.linspace(0, 1, N)

obs_elev = sample_lonlat(args.lon1, args.lat1) + args.h1
tgt_elev = sample_lonlat(args.lon2, args.lat2) + args.h2
profile = [sample_lonlat(lon, lat) for lon, lat in zip(lons, lats)]

visible = True
max_angle = -999
for i in range(1, N):
    dist_ratio = dists[i]
    expected_height = obs_elev + dist_ratio * (tgt_elev - obs_elev)
    dz = profile[i] - expected_height
    if dz > 0:
        visible = False
        break

with open("data/los_result.json", "w") as f:
    json.dump({"visible": visible, "observer": [args.lon1, args.lat1], "target": [args.lon2, args.lat2]}, f)
