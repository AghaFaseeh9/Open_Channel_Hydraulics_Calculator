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

        plt.figure()
        plt.plot(x, bed, label='Bed')
        plt.plot(x, water_surface, label='Water Surface / HGL')
        plt.plot(x, energy_line, label='Energy Line')
        plt.xlabel('Length along channel' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Elevation' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Longitudinal Section')
        plt.legend()
        plt.grid(True)
        plt.show()

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
        plt.figure()
        plt.plot([0, self.b], [0, 0], 'k-')  # Bottom
        plt.plot([0, 0], [0, y], 'b-')  # Left water
        plt.plot([self.b, self.b], [0, y], 'b-')  # Right water
        plt.plot([0, self.b], [y, y], 'c--')  # Water surface
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Rectangular Cross-Section')
        plt.grid(True)
        plt.show()

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
        plt.figure()
        plt.plot([left_x, 0, right_x], [0, y, 0], 'b-')  # Water sides
        plt.plot([left_x, right_x], [0, 0], 'k-')  # Bottom (point)
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Triangular Cross-Section')
        plt.grid(True)
        plt.show()

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
        plt.figure()
        plt.plot([left_x, -self.b/2, self.b/2, right_x], [0, y, y, 0], 'b-')  # Water
        plt.plot([-self.b/2, self.b/2], [y, y], 'c--')  # Water surface
        plt.plot([left_x, right_x], [0, 0], 'k-')  # Bottom
        plt.xlabel('Width' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.ylabel('Depth' + (' (m)' if self.unit == 'si' else ' (ft)'))
        plt.title('Trapezoidal Cross-Section')
        plt.grid(True)
        plt.show()

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
        plt.figure()
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
        plt.show()

# Main function to run the program
def main():
    unit = input("Enter unit system (SI or BG): ").strip()
    shape = input("Enter cross-section shape (rectangular, triangular, trapezoidal, circular): ").strip().lower()
    n = float(input("Enter Manning's n: "))
    S = float(input("Enter bed slope S: "))
    
    if shape == 'rectangular':
        b = float(input("Enter bottom width b: "))
        section = Rectangular(unit, n, S, b)
    elif shape == 'triangular':
        m = float(input("Enter side slope m (horizontal:vertical): "))
        section = Triangular(unit, n, S, m)
    elif shape == 'trapezoidal':
        b = float(input("Enter bottom width b: "))
        m = float(input("Enter side slope m: "))
        section = Trapezoidal(unit, n, S, b, m)
    elif shape == 'circular':
        D = float(input("Enter diameter D: "))
        section = Circular(unit, n, S, D)
    else:
        print("Invalid shape.")
        return

    mode = input("Enter mode (capacity or design): ").strip().lower()
    efficient = False
    if mode == 'design':
        efficient_str = input("Most efficient design? (yes/no): ").strip().lower()
        efficient = efficient_str == 'yes'
        if efficient:
            print("Most efficient configuration:", section.most_efficient())
            # For efficient, we may need to adjust parameters
            # For simplicity, assume user adjusts inputs accordingly, or implement specific
            # Here, we'll proceed with given params, note it's selective

    if mode == 'capacity':
        y = float(input("Enter water depth y: "))
        Q, V, A, P, R = section.compute_Q(y)
    elif mode == 'design':
        Q_target = float(input("Enter target discharge Q: "))
        y_guess = float(input("Enter initial guess for y: "))
        y = section.solve_for_y(Q_target, y_guess)
        Q, V, A, P, R = section.compute_Q(y)
    else:
        print("Invalid mode.")
        return

    # Display inputs and outputs
    print("\nInput Data:")
    print(f"Unit: {unit}")
    print(f"Shape: {shape}")
    print(f"Manning's n: {n}")
    print(f"Slope S: {S}")
    if hasattr(section, 'b'):
        print(f"Bottom width b: {section.b}")
    if hasattr(section, 'm'):
        print(f"Side slope m: {section.m}")
    if hasattr(section, 'D'):
        print(f"Diameter D: {section.D}")
    if mode == 'design':
        print(f"Target Q: {Q_target}")
    print(f"Water depth y: {y}")

    print("\nOutput:")
    print(f"Discharge Q: {Q}" + (' m³/s' if unit.lower() == 'si' else ' cfs'))
    print(f"Velocity V: {V}" + (' m/s' if unit.lower() == 'si' else ' ft/s'))
    print(f"Area A: {A}" + (' m²' if unit.lower() == 'si' else ' ft²'))
    print(f"Wetted Perimeter P: {P}" + (' m' if unit.lower() == 'si' else ' ft'))
    print(f"Hydraulic Radius R: {R}" + (' m' if unit.lower() == 'si' else ' ft'))

    # Plots
    section.plot_cross_section(y)
    section.plot_longitudinal_section(y)

if __name__ == "__main__":
    main()