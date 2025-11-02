import streamlit as st
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import numpy as np

# Function to compute Manning's factor based on unit system
def get_manning_factor(unit):
    if unit.lower() == 'si':
        return 1.0
    elif unit.lower() == 'bg' or unit.lower() == 'fps':
        return 1.486
    else:
        raise ValueError("Invalid unit system. Choose 'SI' or 'BG'.")

# Function to compute g based on unit
def get_g(unit):
    if unit.lower() == 'si':
        return 9.81  # m/s²
    elif unit.lower() == 'bg' or unit.lower() == 'fps':
        return 32.2  # ft/s²
    else:
        raise ValueError("Invalid unit system.")

# Base class for cross-sections
class CrossSection:
    def __init__(self, unit, n, S):
        self.unit = unit
        self.n = n
        self.S = S
        self.factor = get_manning_factor(unit)
        self.g = get_g(unit)

    def compute_Q(self, y):
        A = self.area(y)
        P = self.wetted_perimeter(y)
        R = A / P if P > 0 else 0
        V = (self.factor / self.n) * (R ** (2/3)) * (self.S ** 0.5)
        Q = V * A
        return Q, V, A, P, R

    def solve_for_y(self, Q_target, y_guess=1.0):
        def func(y):
            Q, _, _, _, _ = self.compute_Q(y)
            return Q - Q_target
        y = fsolve(func, y_guess)[0]
        return y

    def plot_cross_section(self, y):
        raise NotImplementedError

    def plot_longitudinal_section(self, y, L=100):
        Q, V, A, _, _ = self.compute_Q(y)
        alpha = V**2 / (2 * self.g)  # Velocity head
        x = np.linspace(0, L, 100)
        bed = -self.S * x
        water_surface = y - self.S * x
        energy_line = y + alpha - self.S * x

        fig = plt.figure()
        plt.plot(x, bed, label='Bed')
        plt.plot(x, water_surface, label='Water Surface / HGL')
        plt.plot(x, energy_line, label='Energy Line')
        plt.xlabel('Length along channel' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Elevation' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Longitudinal Section')
        plt.legend()
        plt.grid(True)
        return fig

# Rectangular cross-section
class Rectangular(CrossSection):
    def __init__(self, unit, n, S, b):
        super().__init__(unit, n, S)
        self.b = b

    def area(self, y):
        return self.b * y

    def wetted_perimeter(self, y):
        return self.b + 2 * y

    def most_efficient(self):
        # For rectangle, most efficient when b = 2y
        return 'b = 2y'

    def plot_cross_section(self, y):
        fig = plt.figure()
        plt.plot([0, self.b], [0, 0], 'k-')  # Bottom
        plt.plot([0, 0], [0, y], 'b-')  # Left water
        plt.plot([self.b, self.b], [0, y], 'b-')  # Right water
        plt.plot([0, self.b], [y, y], 'c--')  # Water surface
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Rectangular Cross-Section')
        plt.grid(True)
        return fig

# Triangular cross-section
class Triangular(CrossSection):
    def __init__(self, unit, n, S, m):  # m = side slope (horizontal:vertical)
        super().__init__(unit, n, S)
        self.m = m

    def area(self, y):
        return self.m * y**2

    def wetted_perimeter(self, y):
        return 2 * y * math.sqrt(1 + self.m**2)

    def most_efficient(self):
        # For triangle, most efficient when m = 1 (45 degrees)
        return 'm = 1 (45 degrees)'

    def plot_cross_section(self, y):
        left_x = -self.m * y
        right_x = self.m * y
        fig = plt.figure()
        plt.plot([left_x, 0, right_x], [0, y, 0], 'b-')  # Water sides
        plt.plot([left_x, right_x], [0, 0], 'k-')  # Bottom (point)
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Triangular Cross-Section')
        plt.grid(True)
        return fig

# Trapezoidal cross-section
class Trapezoidal(CrossSection):
    def __init__(self, unit, n, S, b, m):
        super().__init__(unit, n, S)
        self.b = b
        self.m = m

    def area(self, y):
        return (self.b + self.m * y) * y

    def wetted_perimeter(self, y):
        return self.b + 2 * y * math.sqrt(1 + self.m**2)

    def most_efficient(self):
        # Most efficient trapezoid is half hexagon, m = 1/math.sqrt(3) ≈ 0.577, b = 2*y / math.sqrt(3)
        return 'm = 1/sqrt(3) ≈ 0.577, b = 2y / sqrt(3)'

    def plot_cross_section(self, y):
        left_x = -self.b/2 - self.m * y
        right_x = self.b/2 + self.m * y
        fig = plt.figure()
        plt.plot([left_x, -self.b/2, self.b/2, right_x], [0, y, y, 0], 'b-')  # Water
        plt.plot([-self.b/2, self.b/2], [y, y], 'c--')  # Water surface
        plt.plot([left_x, right_x], [0, 0], 'k-')  # Bottom
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Trapezoidal Cross-Section')
        plt.grid(True)
        return fig

# Circular cross-section
class Circular(CrossSection):
    def __init__(self, unit, n, S, D):
        super().__init__(unit, n, S)
        self.D = D
        self.r = D / 2

    def area(self, y):
        if y >= self.D:
            y = self.D
        if y <= 0:
            return 0
        theta = 2 * math.acos((self.r - y) / self.r)
        return self.r**2 * (theta - math.sin(theta)) / 2

    def wetted_perimeter(self, y):
        if y >= self.D:
            y = self.D
        if y <= 0:
            return 0
        theta = 2 * math.acos((self.r - y) / self.r)
        return self.r * theta

    def most_efficient(self):
        # For max discharge, y ≈ 0.938 D
        return 'y ≈ 0.938 D'

    def plot_cross_section(self, y):
        theta = np.linspace(0, 2*np.pi, 100)
        x = self.r * np.cos(theta)
        z = self.r * np.sin(theta) + self.r  # Shift up
        fig = plt.figure()
        plt.plot(x, z, 'k-')  # Full circle
        water_level = self.r - (self.r - y) if y < self.r else self.r
        plt.axhline(y=water_level, xmin=-self.r, xmax=self.r, color='c', linestyle='--')
        # Approximate water sides
        theta_w = 2 * math.acos((self.r - y) / self.r)
        theta_start = np.pi - theta_w/2
        theta_end = np.pi + theta_w/2
        theta_arc = np.linspace(theta_start, theta_end, 100)
        x_arc = self.r * np.cos(theta_arc)
        z_arc = self.r * np.sin(theta_arc) + self.r
        plt.plot(x_arc, z_arc, 'b-')
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Circular Cross-Section')
        plt.axis('equal')
        plt.grid(True)
        return fig

# Streamlit app
st.title("Open Channel Hydraulics Calculator")

st.markdown("""
This application calculates flow parameters for open channels with various cross-sections using Manning's equation for steady uniform flow.
Select the unit system, shape, and mode to proceed. Ensure all inputs are positive values for accurate results.
""")

# Sidebar for common parameters
with st.sidebar:
    st.header("Common Parameters")
    unit = st.selectbox("Select unit system", ["SI", "BG"], help="SI: Metric units (m, m³/s), BG: British Gravitational (ft, cfs)").strip()
    n = st.number_input("Enter Manning's n", min_value=0.001, value=0.015, step=0.001, help="Roughness coefficient, typically 0.01-0.05 for channels")
    S = st.number_input("Enter bed slope S", min_value=0.0001, value=0.001, step=0.0001, help="Longitudinal slope of the channel bed (dimensionless)")

shape = st.selectbox("Select cross-section shape", ["Rectangular", "Triangular", "Trapezoidal", "Circular"], help="Choose the channel shape").strip().lower()

section = None
unit_length = 'm' if unit.lower() == 'si' else 'ft'
unit_flow = 'm³/s' if unit.lower() == 'si' else 'cfs'
unit_vel = 'm/s' if unit.lower() == 'si' else 'ft/s'
unit_area = 'm²' if unit.lower() == 'si' else 'ft²'

if shape == 'rectangular':
    b = st.number_input(f"Enter bottom width b ({unit_length})", min_value=0.01, value=5.0, step=0.1, help="Width of the rectangular channel bottom")
    if b > 0 and n > 0 and S > 0:
        section = Rectangular(unit, n, S, b)
    else:
        st.warning("Please enter positive values for all parameters.")
elif shape == 'triangular':
    m = st.number_input("Enter side slope m (horizontal:vertical)", min_value=0.1, value=1.0, step=0.1, help="Slope ratio, e.g., 1 for 45 degrees")
    if m > 0 and n > 0 and S > 0:
        section = Triangular(unit, n, S, m)
    else:
        st.warning("Please enter positive values for all parameters.")
elif shape == 'trapezoidal':
    col1, col2 = st.columns(2)
    with col1:
        b = st.number_input(f"Enter bottom width b ({unit_length})", min_value=0.01, value=5.0, step=0.1)
    with col2:
        m = st.number_input("Enter side slope m", min_value=0.1, value=1.5, step=0.1)
    if b > 0 and m > 0 and n > 0 and S > 0:
        section = Trapezoidal(unit, n, S, b, m)
    else:
        st.warning("Please enter positive values for all parameters.")
elif shape == 'circular':
    D = st.number_input(f"Enter diameter D ({unit_length})", min_value=0.01, value=2.0, step=0.1, help="Diameter of the circular channel")
    if D > 0 and n > 0 and S > 0:
        section = Circular(unit, n, S, D)
    else:
        st.warning("Please enter positive values for all parameters.")

if section is not None:
    mode = st.selectbox("Select mode", ["Capacity", "Design"], help="Capacity: Compute Q from y; Design: Compute y from Q").strip().lower()
    efficient = False
    if mode == 'design':
        efficient_str = st.selectbox("Most efficient design?", ["No", "Yes"], help="If yes, displays guidelines for most hydraulic efficient section").strip().lower()
        efficient = efficient_str == 'yes'
        if efficient:
            st.info(f"Most efficient configuration for {shape.capitalize()}: {section.most_efficient()} \nAdjust inputs accordingly and recalculate.")

    y = None
    Q_target = None
    if mode == 'capacity':
        y = st.number_input(f"Enter water depth y ({unit_length})", min_value=0.01, value=1.0, step=0.1, help="Flow depth in the channel")
    elif mode == 'design':
        Q_target = st.number_input(f"Enter target discharge Q ({unit_flow})", min_value=0.01, value=10.0, step=1.0, help="Desired flow rate")
        y_guess = st.number_input(f"Enter initial guess for y ({unit_length})", min_value=0.01, value=1.0, step=0.1, help="Starting guess for solver")
        if st.button("Calculate", help="Click to compute the flow depth"):
            with st.spinner("Calculating..."):
                try:
                    y = section.solve_for_y(Q_target, y_guess)
                    if y <= 0:
                        st.error("Solver returned non-positive depth. Try a different initial guess.")
                        y = None
                except Exception as e:
                    st.error(f"Error in calculation: {str(e)}")
                    y = None

    if y is not None and y > 0:
        Q, V, A, P, R = section.compute_Q(y)

        # Display inputs and outputs in expanders
        with st.expander("Input Data", expanded=True):
            st.write(f"**Unit System:** {unit.upper()}")
            st.write(f"**Shape:** {shape.capitalize()}")
            st.write(f"**Manning's n:** {n:.3f}")
            st.write(f"**Slope S:** {S:.4f}")
            if hasattr(section, 'b'):
                st.write(f"**Bottom width b ({unit_length}):** {section.b:.2f}")
            if hasattr(section, 'm'):
                st.write(f"**Side slope m:** {section.m:.2f}")
            if hasattr(section, 'D'):
                st.write(f"**Diameter D ({unit_length}):** {section.D:.2f}")
            if mode == 'design':
                st.write(f"**Target Q ({unit_flow}):** {Q_target:.2f}")
            st.write(f"**Water depth y ({unit_length}):** {y:.2f}")

        with st.expander("Output Results", expanded=True):
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Discharge Q", f"{Q:.2f} {unit_flow}")
                st.metric("Velocity V", f"{V:.2f} {unit_vel}")
            with col2:
                st.metric("Area A", f"{A:.2f} {unit_area}")
                st.metric("Wetted Perimeter P", f"{P:.2f} {unit_length}")
            with col3:
                st.metric("Hydraulic Radius R", f"{R:.2f} {unit_length}")

        # Plots in tabs
        tab1, tab2 = st.tabs(["Cross-Section Plot", "Longitudinal Section Plot"])
        with tab1:
            cross_fig = section.plot_cross_section(y)
            st.pyplot(cross_fig)
        with tab2:
            long_fig = section.plot_longitudinal_section(y)
            st.pyplot(long_fig)
else:
    if shape:
        st.info("Please enter valid positive parameters for the selected shape to proceed.")