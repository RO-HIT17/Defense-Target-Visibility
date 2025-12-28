# Defense Visibility & Line of Sight System

A geospatial web application for defense planning, allowing users to simulate viewshed (visibility) and line-of-sight (LOS) analysis over terrain using digital elevation models (DEMs).

## Features
- **Viewshed Simulation:** Compute and visualize areas visible from an observer point at a given height and location.
- **Line of Sight Checker:** Check if a direct line of sight exists between two points (observer and target) with specified heights.
- **Interactive Maps:** Results are displayed on interactive Folium maps in the browser.
- **Backend Scripts:** FastAPI backend for geospatial computation; Python scripts for raster viewshed analysis.

## Project Structure
```
DefenseVisiblity/
├── defense_los/
│   ├── app.py                # Streamlit frontend app
│   ├── compute_viewshed.py   # Standalone viewshed computation script
│   ├── ...                   # Other backend files
├── data/                     # DEMs and output files
├── README.md                 # Project documentation
```

## Setup
1. **Install dependencies:**
   - Python 3.8+
   - Install required packages:
     ```bash
     pip install streamlit folium streamlit-folium requests rasterio numpy scipy fastapi uvicorn
     ```
2. **Prepare DEM data:**
   - Place your GeoTIFF DEM file in the `data/` folder.

## Running the App
1. **Start the backend (FastAPI):**
   - (Assuming you have a FastAPI backend exposing `/viewshed` and `/los` endpoints)
   - Example:
     ```bash
     uvicorn backend:app --reload --port 8000
     ```
2. **Start the frontend (Streamlit):**
   ```bash
   streamlit run defense_los/app.py
   ```
3. **Open your browser:**
   - Visit [http://localhost:8501](http://localhost:8501)

## Standalone Viewshed Script
- Run viewshed analysis directly from CLI:
  ```bash
  python defense_los/compute_viewshed.py --lon <lon> --lat <lat> --height <h> --maxdist <m> --dem data/demo_aoi.tif
  ```
- Output is a GeoJSON file with visible polygons.

## Notes
- For best results, use high-resolution DEMs.
- The app caches results for faster repeated queries.
- Morphological smoothing in viewshed script requires `scipy`.
