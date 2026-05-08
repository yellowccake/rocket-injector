import math

# =====================
# INPUT VALUES
# =====================

F = 40_000                 # thrust [N]
g0 = 9.81                  # gravity [m/s^2]
pc = 80e5                  # chamber pressure [Pa]

ROF = 3.5                  # oxidizer/fuel ratio

c_ideal = 3084.1           # ideal exhaust velocity from CEA [m/s]
cstar_ideal = 1845.2       # ideal c* from CEA [m/s]

eta_isp = 0.88             # impulse efficiency
eta_cf = 0.96              # nozzle / CF efficiency

contraction_ratio = 2.5    # A_face / A_throat


# =====================
# CALCULATIONS
# =====================

# Real specific impulse
Isp_ideal = c_ideal / g0
Isp_real = Isp_ideal * eta_isp

# Total mass flow
mdot_total = F / (Isp_real * g0)

# Fuel and oxidizer mass flow
mdot_fuel = mdot_total / (ROF + 1)
mdot_oxidizer = ROF * mdot_fuel

# c* efficiency
eta_cstar = eta_isp / eta_cf

# Real c*
cstar_real = cstar_ideal * eta_cstar

# Throat area
A_throat = cstar_real * mdot_total / pc

# Throat diameter
d_throat = 2 * math.sqrt(A_throat / math.pi)

# Face plate / injector face area
A_face = contraction_ratio * A_throat

# Face plate diameter
d_face = 2 * math.sqrt(A_face / math.pi)


# =====================
# OUTPUT
# =====================

print("===== INPUTS =====")
print(f"Thrust = {F:.1f} N")
print(f"Pc = {pc/1e5:.1f} bar")
print(f"ROF = {ROF:.2f}")
print(f"c ideal = {c_ideal:.2f} m/s")
print(f"c* ideal = {cstar_ideal:.2f} m/s")

print("\n===== RESULTS =====")
print(f"Isp ideal = {Isp_ideal:.4f} s")
print(f"Isp real = {Isp_real:.4f} s")

print(f"Total mass flow = {mdot_total:.4f} kg/s")
print(f"Fuel mass flow = {mdot_fuel:.4f} kg/s")
print(f"Oxidizer mass flow = {mdot_oxidizer:.4f} kg/s")

print(f"eta c* = {eta_cstar:.4f}")
print(f"c* real = {cstar_real:.4f} m/s")

print(f"A throat = {A_throat:.6e} m^2")
print(f"A throat = {A_throat*1e6:.4f} mm^2")

print(f"Throat diameter = {d_throat*1000:.4f} mm")

print(f"Face plate area = {A_face*1e6:.4f} mm^2")
print(f"Face plate diameter = {d_face*1000:.4f} mm")