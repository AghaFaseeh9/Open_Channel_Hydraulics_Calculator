# 🌊 Open Channel Hydraulics Calculator

A **Streamlit web application** that computes and visualizes **open channel flow parameters** for different cross-sectional shapes using **Manning’s Equation** for steady uniform flow.  
This tool is designed for **civil and hydraulic engineering students** to analyze flow behavior in rectangular, triangular, trapezoidal, and circular channels with both **SI** and **British Gravitational (BG/FPS)** unit systems.

---

## 🚀 Features

✅ **Supports multiple cross-sections**
- Rectangular  
- Triangular  
- Trapezoidal  
- Circular  

✅ **Two calculation modes**
- **Capacity Mode:** Computes discharge (Q) and velocity (V) for a given flow depth (y).  
- **Design Mode:** Calculates required flow depth (y) for a target discharge (Q).  

✅ **Unit System Selection**
- Choose between **SI (m, m³/s)** and **BG (ft, cfs)** systems.

✅ **Graphical Visualizations**
- **Cross-section plot:** Displays water depth and geometry.  
- **Longitudinal section plot:** Shows bed profile, hydraulic grade line, and energy line.

✅ **Interactive Inputs**
- Built with **Streamlit widgets** for real-time calculations and plotting.  
- Includes error handling and user guidance for realistic parameters.

✅ **Most Efficient Channel Option**
- Displays conditions for hydraulically efficient section geometry.

---

## 🧮 Formula Used: Manning’s Equation

**Manning’s Equation:**  

Q = (1 / n) × A × R^(2/3) × S^(1/2)




Where:  
- \( Q \) = Discharge  
- \( n \) = Manning’s roughness coefficient  
- \( A \) = Flow area  
- \( R \) = Hydraulic radius (A/P)  
- \( S \) = Channel bed slope  

For BG (FPS) units, the equation includes a factor of **1.486**.

---
