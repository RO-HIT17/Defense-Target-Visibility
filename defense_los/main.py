from fastapi import FastAPI, Query
from fastapi.responses import JSONResponse
import subprocess
import json
import os
import uvicorn

app = FastAPI(title="Defense Visibility API")

@app.get("/viewshed")
def compute_viewshed(
    lon: float = Query(..., ge=70.99986111111112, le=72.0001388888889),
    lat: float = Query(..., ge=33.99986111111111, le=35.000138888888884),
    height: float = Query(50.0),
    maxdist: int = Query(5000),
):
    """Compute viewshed for given observer coordinates and height."""
    cmd = [
        "python", "compute_viewshed.py",
        "--lon", str(lon),
        "--lat", str(lat),
        "--height", str(height),
        "--maxdist", str(maxdist),
        "--dem", "data/demo_aoi.tif"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return JSONResponse(status_code=500, content={"error": result.stderr})

    with open("data/viewshed_py.geojson") as f:
        geojson = json.load(f)
    return geojson


@app.get("/los")
def compute_los(
    lon1: float = Query(...),
    lat1: float = Query(...),
    lon2: float = Query(...),
    lat2: float = Query(...),
    h1: float = Query(30.0),
    h2: float = Query(10.0),
):
    """Compute line-of-sight between two points."""
    cmd = [
        "python", "compute_los.py",
        "--lon1", str(lon1), "--lat1", str(lat1),
        "--lon2", str(lon2), "--lat2", str(lat2),
        "--h1", str(h1), "--h2", str(h2),
        "--dem", "data/demo_aoi.tif"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return JSONResponse(status_code=500, content={"error": result.stderr})

    with open("data/los_result.json") as f:
        data = json.load(f)
    return data


if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
