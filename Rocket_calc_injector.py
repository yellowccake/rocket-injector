import math
import cantera as ct

# INPUT VALUES
# This is another test (Diego)
# Test From JP
#stop testing guys, it works lol

F = 40_000                 # thrust [N]
g0 = 9.81                  # gravity [m/s^2]
pc = 80e5                  # chamber pressure [Pa]



ROF = 3.25                 # oxidizer/fuel ratio

c_ideal = 3100.3           # ideal exhaust velocity from CEA [m/s]
cstar_ideal = 1866       # ideal c* from CEA [m/s]

eta_isp = 0.88             # impulse efficiency
eta_cf = 0.96              # nozzle / CF efficiency

contraction_ratio = 2.5    # A_face / A_throat (Source: Alex)

T_fuel = 250               # CH4 inlet temperature [K] (fuel=CH4)
T_Lox = 96                  # LOX inlet temperature [K]

# PHASE 1: Engine performance & mass flow
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

# PHASE 2: Throat & Injector face sizing
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



N_elements = 19
n_inlets = 4

delta_p_inj = 16e5

rho_LOX = 1157.9

alpha_half_deg = 50
alpha_full_deg = 2 * alpha_half_deg

alpha_half_rad = math.radians(alpha_half_deg)

phi = 2 / (math.tan(alpha_half_rad)**2 + 2)

A_bazarov = math.sqrt(
    2 * (1 - phi)**2 / phi**3
)

mu_swirl = phi * math.sqrt(
    phi / (2 - phi)
)

mdot_LOX_el = mdot_oxidizer / N_elements

Rn_LOX = 0.475 * math.sqrt(
    mdot_LOX_el /
    (
        mu_swirl *
        math.sqrt(rho_LOX * delta_p_inj)
    )
)

Dn_LOX = 2 * Rn_LOX

Rin = 2 * Rn_LOX

rin = math.sqrt(
    (Rin * Rn_LOX) /
    (n_inlets * A_bazarov)
)

din = 2 * rin

l_in = 5 * rin

l_n = 1 * Rn_LOX

R_s = Rin + rin

D_s = 2 * R_s

l_s = 2 * Rin

nu_LOX = 1.8160e-7
xi_in = 0.77

Re_in = 0.637 * mdot_LOX_el / (
    math.sqrt(n_inlets) * rin * rho_LOX * nu_LOX
)

lambda_f = 0.3164 / (Re_in ** 0.25)

A_eq = (Rin * Rn_LOX) / (
    n_inlets * rin**2 + (lambda_f / 2) * Rin * (Rin - Rn_LOX)
)

xi_total = xi_in + lambda_f * (R_s / Rin)

mu_i = mu_swirl / math.sqrt(
    1 + xi_total * mu_swirl**2 * A_eq**2 / ((Rin / Rn_LOX)**2)
)

def calc_phi_from_A(A):
    low = 1e-6
    high = 1 - 1e-6

    for _ in range(50):
        mid = (low + high) / 2
        value = A**2 * mid**3 - 2 * (1 - mid)**2

        if value > 0:
            high = mid
        else:
            low = mid

    return (low + high) / 2

tolerance = 1e-4
max_iterations = 20


A_current = A_eq
mu_i_current = mu_i
Rn_LOX_corrected = Rn_LOX

for iteration in range(max_iterations):
    Rn_LOX_corrected = 0.475 * math.sqrt(
        mdot_LOX_el / (mu_i_current * math.sqrt(rho_LOX * delta_p_inj))
    )

    A_new = (Rin * Rn_LOX_corrected) / (n_inlets * rin**2)

    phi_new = calc_phi_from_A(A_new)

    mu_eq_new = phi_new * math.sqrt(phi_new / (2 - phi_new))

    mu_i_new = mu_eq_new / math.sqrt(
        1 + xi_total * mu_eq_new**2 * A_new**2 / ((Rin / Rn_LOX_corrected)**2)
    )

    if abs(A_new - A_current) < tolerance:
        break

    A_current = A_new
    mu_i_current = mu_i_new

A_corrected = A_current
mu_i = mu_i_current
Dn_LOX_corrected = 2 * Rn_LOX_corrected

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
print("\n===== PHASE 3 : LOX SWIRL INJECTOR SIZING =====")


print(f"Injector elements = {N_elements}")
print(f"Tangential inlets per element = {n_inlets}")
print(f"LOX mass flow per element = {mdot_LOX_el:.4f} kg/s")

print(f"Spray half-angle = {alpha_half_deg:.1f} deg")
print(f"Spray full cone angle = {alpha_full_deg:.1f} deg")

print("\n--- Ideal Bazarov values ---")
print(f"A_bazarov ideal = {A_bazarov:.4f} [-]")
print(f"phi = {phi:.4f} [-]")
print(f"mu_swirl ideal = {mu_swirl:.4f} [-]")
print(f"LOX nozzle radius ideal Rn = {Rn_LOX*1000:.3f} mm")
print(f"LOX nozzle diameter ideal Dn = {Dn_LOX*1000:.3f} mm")

print("\n--- Real correction values ---")
print(f"Re_in = {Re_in:.1f} [-]")
print(f"lambda_f = {lambda_f:.5f} [-]")
print(f"A_eq = {A_eq:.4f} [-]")
print(f"xi_total = {xi_total:.4f} [-]")
print(f"mu_i corrected = {mu_i:.4f} [-]")
print(f"LOX nozzle radius corrected Rn = {Rn_LOX_corrected*1000:.3f} mm")
print(f"LOX nozzle diameter corrected Dn = {Dn_LOX_corrected*1000:.3f} mm")

print(f"Iterations = {iteration + 1}")
print(f"A corrected = {A_corrected:.4f} [-]")
print(f"mu_i corrected = {mu_i:.4f} [-]")
print(f"LOX nozzle diameter corrected Dn = {Dn_LOX_corrected*1000:.3f} mm")