
import math
import numpy as np
import rasterio
from rasterio.features import shapes
import json
from shapely.geometry import shape, mapping
import fiona
from fiona.crs import from_epsg
from rasterio.transform import from_origin


DEM_PATH = r"data\demo_crop.tif"   
OUT_TIF = r"data\viewshed_py.tif"
OUT_GEOJSON = r"data\viewshed_py.geojson"


OBSERVER_LON = 71.385
OBSERVER_LAT = 34.189

OBSERVER_HEIGHT_M = 50.0      
MAX_DISTANCE_M = 5000        
AZIMUTH_STEPS = 720          
STEP_M = None                



METERS_PER_DEG_LAT = 111132.92  
def meters_per_deg_lon(lat_deg):
    return 111412.84 * math.cos(math.radians(lat_deg)) - 93.5 * math.cos(3*math.radians(lat_deg))

def run():
    
    with rasterio.open(DEM_PATH) as src:
        dem = src.read(1, masked=True)   
        transform = src.transform
        crs = src.crs
        width = src.width
        height = src.height

    
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
        raise SystemExit("Observer point falls on NoData in DEM. Pick a nearby valid coordinate.")

    
    visible = np.zeros_like(dem, dtype=np.uint8)  

    
    
    mpd_lat = METERS_PER_DEG_LAT
    mpd_lon = meters_per_deg_lon(OBSERVER_LAT)
    
    pix_size_x = transform.a   
    pix_size_y = -transform.e  
    meters_per_pixel = math.hypot(pix_size_x * mpd_lon, pix_size_y * mpd_lat)
    if STEP_M is None:
        step_m = meters_per_pixel * 0.9
    else:
        step_m = STEP_M

    
    max_steps = int(math.ceil(MAX_DISTANCE_M / step_m))
    
    two_pi = 2 * math.pi
    for i, theta in enumerate(np.linspace(0, two_pi, AZIMUTH_STEPS, endpoint=False)):
        
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

    print("Done.")
    print("Outputs:")
    print(" - TIFF:", OUT_TIF)
    print(" - GeoJSON:", OUT_GEOJSON)

if __name__ == "__main__":
    run()
