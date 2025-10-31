import streamlit as st
import folium
from streamlit_folium import st_folium
import requests

API_URL = "http://localhost:8000"


@st.cache_data(show_spinner=False)
def get_viewshed(lon: float, lat: float, height: float, maxdist: int):
    """Fetch viewshed geojson from backend and cache by input parameters."""
    res = requests.get(f"{API_URL}/viewshed", params={
        "lon": lon, "lat": lat, "height": height, "maxdist": maxdist
    }, timeout=15)
    if res.status_code != 200:
        raise RuntimeError(res.text)
    return res.json()

@st.cache_data(show_spinner=False)
def get_los(lon1: float, lat1: float, lon2: float, lat2: float, h1: float, h2: float):
    """Fetch LOS result from backend and cache by input parameters."""
    res = requests.get(f"{API_URL}/los", params={
        "lon1": lon1, "lat1": lat1, "lon2": lon2, "lat2": lat2, "h1": h1, "h2": h2
    }, timeout=15)
    if res.status_code != 200:
        raise RuntimeError(res.text)
    return res.json()

st.set_page_config(page_title="Defense Visibility System", layout="wide")
st.title("🛰️ Defense Visibility & Line of Sight System")

tab1, tab2 = st.tabs(["Defense Visibility", "Line of Sight"])


for key in ["viewshed_geojson", "viewshed_center", "los_params", "los_visible"]:
    if key not in st.session_state:
        st.session_state[key] = None





with tab1:
    st.subheader("Defense Visibility Simulation")

    with st.form("viewshed_form"):
        lon = st.number_input("Observer Longitude", 70.9998, 72.0001, 71.985, 0.0001)
        lat = st.number_input("Observer Latitude", 33.9998, 35.0001, 34.189, 0.0001)
        height = st.slider("Observer Height (m)", 0, 100, 40)
        maxdist = st.slider("Max Distance (m)", 1000, 10000, 5000, step=500)
        submitted = st.form_submit_button("Compute Viewshed")

    if submitted:
        try:
            st.info("Fetching data from backend...")
            geojson = get_viewshed(lon, lat, height, maxdist)
            st.session_state.viewshed_geojson = geojson
            st.session_state.viewshed_center = (lat, lon)
            st.success("✅ Viewshed generated successfully!")
        except Exception as e:
            st.error(f"Error connecting to backend: {e}")

    
    if st.session_state.viewshed_geojson and st.session_state.viewshed_center:
        lat_c, lon_c = st.session_state.viewshed_center
        m = folium.Map(location=[lat_c, lon_c], zoom_start=11)
        folium.Marker([lat_c, lon_c], tooltip="Observer", icon=folium.Icon(color="red")).add_to(m)
        folium.GeoJson(
            st.session_state.viewshed_geojson,
            style_function=lambda x: {"fillColor": "green", "color": "green", "fillOpacity": 0.4}
        ).add_to(m)
        st_folium(m, width=1100, height=600, key="viewshed_map")





with tab2:
    st.subheader("Line of Sight Checker")

    with st.form("los_form"):
        c1, c2 = st.columns(2)
        with c1:
            lon1 = st.number_input("Observer Longitude", 70.9998, 72.0001, 71.985, 0.0001)
            lat1 = st.number_input("Observer Latitude", 33.9998, 35.0001, 34.189, 0.0001)
            h1 = st.slider("Observer Height (m)", 0, 100, 40)
        with c2:
            lon2 = st.number_input("Target Longitude", 70.9998, 72.0001, 71.983, 0.0001)
            lat2 = st.number_input("Target Latitude", 33.9998, 35.0001, 34.182, 0.0001)
            h2 = st.slider("Target Height (m)", 0, 100, 20)
        submitted = st.form_submit_button("Check Line of Sight")

    if submitted:
        try:
            st.info("Checking LOS with backend...")
            data = get_los(lon1, lat1, lon2, lat2, h1, h2)
            visible = data.get("visible", False)
            st.session_state.los_params = {"lon1": lon1, "lat1": lat1, "lon2": lon2, "lat2": lat2}
            st.session_state.los_visible = visible
            st.success("✅ LOS data received!")
        except Exception as e:
            st.error(f"Error connecting to backend: {e}")

    
    if st.session_state.los_params is not None:
        p = st.session_state.los_params
        color = "green" if st.session_state.los_visible else "red"
        m = folium.Map(location=[(p["lat1"] + p["lat2"]) / 2, (p["lon1"] + p["lon2"]) / 2], zoom_start=12)
        folium.Marker([p["lat1"], p["lon1"]], tooltip="Observer", icon=folium.Icon(color="blue")).add_to(m)
        folium.Marker([p["lat2"], p["lon2"]], tooltip="Target", icon=folium.Icon(color="orange")).add_to(m)
        folium.PolyLine([[p["lat1"], p["lon1"]], [p["lat2"], p["lon2"]]], color=color, weight=4).add_to(m)
        st_folium(m, width=1100, height=600, key="los_map")
        if st.session_state.los_visible is not None:
            st.success(f"Line of Sight: {'VISIBLE ✅' if st.session_state.los_visible else 'BLOCKED ❌'}")
