import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(layout="wide", page_title="Defense 3D Visibility")

st.title("3D Defense Visibility Viewer (CesiumJS)")

with open("templates/cesium_map.html", "r", encoding="utf-8") as f:
    html_code = f.read()

# Render the Cesium viewer
components.html(html_code, height=800, scrolling=False)
