import math
import cantera as ct

# INPUT VALUES
# This is another test (Diego)
# Test From JP

F = 40_000                 # thrust [N]
g0 = 9.81                  # gravity [m/s^2]
pc = 80e5                  # chamber pressure [Pa]



ROF = 3.25                 # oxidizer/fuel ratio

c_ideal = 3100.3           # ideal exhaust velocity from CEA [m/s]
cstar_ideal =1866       # ideal c* from CEA [m/s]

eta_isp = 0.88             # impulse efficiency
eta_cf = 0.96              # nozzle / CF efficiency

contraction_ratio = 2.5    # A_face / A_throat (Source: Alex)

T_fuel = 250               # CH4 inlet temperature [K] (fuel=CH4)
T_Lox = 96                  # LOX inlet temperature [K]


# CANTERA CHECK


gas = ct.Solution("gri30.yaml")

MW_CH4 = 16.043            # kg/kmol
MW_O2 = 31.999             # kg/kmol


CH4_moles = 1.0
O2_moles = ROF * MW_CH4 / MW_O2

# Simple Cantera approximation:
# one mixed reactant temperature.
# This is NOT exactly the same as NASA CEA with separate CH4/LOX temperatures.
T_mix = 300

gas.TPX = T_mix, pc, f"CH4:{CH4_moles}, O2:{O2_moles}"
gas.equilibrate("HP")

T_chamber_cantera = gas.T
gamma_cantera = gas.cp_mass / gas.cv_mass
MW_cantera = gas.mean_molecular_weight



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


# OUTPUT


print("===== INPUTS =====")
print(f"Thrust = {F:.1f} N")
print(f"Pc = {pc/1e5:.1f} bar")
print(f"ROF = {ROF:.2f}")
print(f"c ideal = {c_ideal:.2f} m/s")
print(f"c* ideal = {cstar_ideal:.2f} m/s")
print(f"Fuel temp = {T_fuel:.1f} K")
print(f"Ox temp = {T_Lox:.1f} K")

print("\n===== CANTERA CHECK =====")
print(f"Cantera CH4 moles = {CH4_moles:.4f}")
print(f"Cantera O2 moles = {O2_moles:.4f}")
print(f"Cantera chamber temperature = {T_chamber_cantera:.2f} K")
print(f"Cantera gamma = {gamma_cantera:.4f}")
print(f"Cantera molecular weight = {MW_cantera:.4f} kg/kmol")

print("\n===== HAND / CEA-BASED RESULTS =====")
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