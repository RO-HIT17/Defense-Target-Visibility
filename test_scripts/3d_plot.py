import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
import pyvista as pv

DEM_PATH = "data/demo_crop.tif"
VIEW_PATH = "data/viewshed_py.tif"

# Load DEM
with rasterio.open(DEM_PATH) as src:
    dem = src.read(1)
    dem_transform = src.transform
    dem_profile = src.profile

# Load Viewshed and resample to DEM grid
with rasterio.open(VIEW_PATH) as vs_src:
    viewshed = vs_src.read(1)
    viewshed_resampled = np.empty_like(dem)

    reproject(
        source=viewshed,
        destination=viewshed_resampled,
        src_transform=vs_src.transform,
        src_crs=vs_src.crs,
        dst_transform=dem_transform,
        dst_crs=dem_profile['crs'],
        resampling=Resampling.nearest
    )

# Now both have same shape
print("DEM shape:", dem.shape)
print("Viewshed shape:", viewshed_resampled.shape)

# Coordinates
rows, cols = np.indices(dem.shape)
xs = dem_transform.c + cols * dem_transform.a + rows * dem_transform.b
ys = dem_transform.f + cols * dem_transform.d + rows * dem_transform.e

grid = pv.StructuredGrid(xs, ys, dem)

# Colors
visible_colors = np.zeros((*dem.shape, 3))
visible_colors[viewshed_resampled == 1] = [0, 1, 0]
visible_colors[viewshed_resampled == 0] = [1, 0, 0]

# 3D plot
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars=dem, cmap="terrain", opacity=0.8)
plotter.add_axes()
plotter.add_mesh(pv.Sphere(radius=100, center=(71.985, 34.189, np.nanmax(dem)+100)), color="yellow")
plotter.show()
