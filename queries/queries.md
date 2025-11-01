CREATE DATABASE defense_visibility;
\c defense_visibility;
CREATE EXTENSION postgis;
CREATE EXTENSION postgis_raster;


raster2pgsql -s 4326 -I -C -M data/demo_aoi.tif public.dem | psql -d defense_visibility

-- Observer info
WITH observer AS (
  SELECT
    ST_SetSRID(ST_MakePoint(71.385, 34.189), 4326) AS geom
)
SELECT
  ST_Viewshed(
    rast,                          -- DEM raster
    1,                             -- Band number
    geom,                          -- Observer point
    50.0,                          -- Observer height (m)
    0.0,                           -- Target height (optional)
    5000.0,                        -- Max distance in meters
    0.0                            -- Earth curvature correction (0 = off)
  ) AS viewshed_raster
INTO viewshed_result
FROM dem, observer
LIMIT 1;


gdal_translate "PG:dbname=defense_visibility schema=public table=viewshed_result mode=2" data/viewshed_postgis.tif


CREATE TABLE viewshed_poly AS
SELECT
  (ST_DumpAsPolygons(rast)).*
FROM viewshed_result;

-- Export polygons with visible=1
COPY (
  SELECT ST_AsGeoJSON(geom) AS geometry
  FROM viewshed_poly
  WHERE val = 1
) TO 'C:\path\data\viewshed_postgis.geojson';
