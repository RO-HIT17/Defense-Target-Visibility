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
            lon1 = st.number_input("Observer Longitude", 70.9998, 72.0001, 71.985, 0.0001, key="los_lon1")
            lat1 = st.number_input("Observer Latitude", 33.9998, 35.0001, 34.189, 0.0001, key="los_lat1")
            h1 = st.slider("Observer Height (m)", 0, 100, 40, key="los_h1")
        with c2:
            lon2 = st.number_input("Target Longitude", 70.9998, 72.0001, 71.983, 0.0001, key="los_lon2")
            lat2 = st.number_input("Target Latitude", 33.9998, 35.0001, 34.182, 0.0001, key="los_lat2")
            h2 = st.slider("Target Height (m)", 0, 100, 20, key="los_h2")
        submitted = st.form_submit_button("Check Line of Sight")

    if submitted:
        try:
            st.info("Checking LOS with backend...")
            data = get_los(lon1, lat1, lon2, lat2, h1, h2)
            st.session_state.los_params = data
            st.success("✅ LOS data received!")
        except Exception as e:
            st.error(f"Error connecting to backend: {e}")

    if st.session_state.los_params:
        p = st.session_state.los_params
        visible = p.get("visible", False)
        blocking_point = p.get("blocking_point")
        obs = p.get("observer", {})
        tgt = p.get("target", {})

        color = "green" if visible else "red"
        m = folium.Map(location=[(obs["lat"] + tgt["lat"]) / 2, (obs["lon"] + tgt["lon"]) / 2],zoom_start=12        )

        # Observer Marker
        folium.Marker(
            [obs["lat"], obs["lon"]],
            tooltip=(
                f"Observer\n"
                f"Ground: {obs['ground']:.1f} m\n"
                f"Height: {obs['height']:.1f} m\n"
                f"Total Elev: {obs['elev']:.1f} m"
            ),
            icon=folium.Icon(color="blue", icon="user")
        ).add_to(m)

        # Target Marker
        folium.Marker(
            [tgt["lat"], tgt["lon"]],
            tooltip=(
                f"Target\n"
                f"Ground: {tgt['ground']:.1f} m\n"
                f"Height: {tgt['height']:.1f} m\n"
                f"Total Elev: {tgt['elev']:.1f} m"
            ),
            icon=folium.Icon(color="orange", icon="flag")
        ).add_to(m)

        # Line between
        folium.PolyLine(
            [[obs["lat"], obs["lon"]], [tgt["lat"], tgt["lon"]]],
            color=color, weight=4
        ).add_to(m)

        # Blocking point marker (if any)
        if blocking_point:
            folium.Marker(
                [blocking_point["lat"], blocking_point["lon"]],
                tooltip=(
                    f"Blocking Point\n"
                    f"Elev: {blocking_point['elev']:.1f} m"
                ),
                icon=folium.Icon(color="red", icon="remove-sign")
            ).add_to(m)

        st_folium(m, width=1100, height=600, key="los_map")

        if visible:
            st.success("✅ Line of Sight: CLEAR — Target is visible")
        else:
            st.error("❌ Line of Sight: BLOCKED — Obstruction detected!")
            if blocking_point:
                st.info(
                    f"Blocking Point: ({blocking_point['lat']:.5f}, {blocking_point['lon']:.5f}), "
                    f"Elevation: {blocking_point['elev']:.1f} m"
                )

        st.markdown("---")
        st.markdown(
            f"""
            **Observer Elevation:** {obs['elev']:.1f} m  
            **Target Elevation:** {tgt['elev']:.1f} m  
            **Status:** {"✅ Visible" if visible else "❌ Blocked"}
            """
        )
