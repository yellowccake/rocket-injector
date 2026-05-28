"""
SWIRL INJECTOR DESIGN  –  Reproduction of Figure 32
Reference: Yang (2004), "Liquid Rocket Thrust Chambers"
           Chapter: Bazarov, Yang & Puri,
           "Design and Dynamics of Jet and Swirl Injectors"

Figure 32 shows the effect of the geometric characteristic
parameter A on six key dimensionless injector parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ─────────────────────────────────────────────────────────
#  NOMENCLATURE
#   A      : geometric characteristic parameter       Eq. (48)
#   phi    : coefficient of passage fullness (φ)      Eq. (46)
#   mu     : mass flow (discharge) coefficient (μ)    Eq. (62)
#   U_un   : non-dim. azimuthal velocity in nozzle    Eq. (72)
#   U_an   : non-dim. axial velocity in nozzle        Eq. (73)
#   alpha_n: spray cone half-angle at nozzle exit     Eq. (74)
#   h_bar  : non-dim. vortex-core radius at VC head   Eq. (75)
#            h_bar = U_ue = sqrt(a),  a = r_mk_bar^2
# ─────────────────────────────────────────────────────────


# ── 1.  Helper: solve for phi given A ──────────────────────────────────────
def calc_phi(A):
    """
    Solve for phi from the combined Eqs. (61) and (62):
        A^2 * phi^3 = 2*(1 - phi)^2
    """
    fun = lambda p: A**2 * p**3 - 2 * (1 - p)**2
    return brentq(fun, 1e-6, 1 - 1e-6)


# ── 2.  Sweep over A and solve for phi ─────────────────────────────────────
N      = 1000
A_vec  = np.linspace(0.05, 14, N)
phi    = np.array([calc_phi(A) for A in A_vec])


# ── 3.  Compute all parameters from phi ────────────────────────────────────

# Eq. (62)  – mass flow / discharge coefficient
mu    = phi * np.sqrt(phi / (2 - phi))

# Eq. (72)  – non-dim. azimuthal velocity in nozzle
U_un  = np.sqrt(2 * (1 - phi) / (2 - phi))

# Eq. (73)  – non-dim. axial velocity in nozzle
U_an  = np.sqrt(phi / (2 - phi))

# Eq. (75)  – h_bar = r_mk_bar = sqrt(a),  a = 2(1-phi)^2/(2-phi)
a_par = 2 * (1 - phi)**2 / (2 - phi)
h_bar = np.sqrt(a_par)

# Eq. (74)  – spray cone half-angle at nozzle exit [degrees]
alpha_n = np.degrees(np.arctan(np.sqrt(2 * (1 - phi) / phi)))


# ── 4.  Plot – reproduce Fig. 32 ───────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(10, 6.5), facecolor='white')

# ---- Left axis : dimensionless parameters ----
h1, = ax1.plot(A_vec, U_un,  'b-',  linewidth=2.0, label=r'$\bar{U}_{un}$  azimuthal vel. (nozzle)')
h2, = ax1.plot(A_vec, h_bar, 'k-',  linewidth=2.0, label=r'$\bar{h} = \sqrt{a}$  vortex-core radius (VC head)')
h3, = ax1.plot(A_vec, U_an,  'r-',  linewidth=2.0, label=r'$\bar{U}_{an}$  axial vel. (nozzle)')
h4, = ax1.plot(A_vec, phi,   'm--', linewidth=1.8, label=r'$\phi$  passage fullness coefficient')
h5, = ax1.plot(A_vec, mu,    'g--', linewidth=1.8, label=r'$\mu$  mass flow (discharge) coefficient')

ax1.set_xlabel('Geometric Characteristic Parameter,  A', fontsize=13, fontweight='bold')
ax1.set_ylabel('Dimensionless Parameters',               fontsize=13, fontweight='bold')
ax1.set_xlim([0, 14])
ax1.set_ylim([0, 1.2])
ax1.set_xticks(np.arange(0, 15, 2))
ax1.set_yticks(np.arange(0, 1.3, 0.2))
ax1.tick_params(axis='both', labelsize=11)
ax1.grid(True, alpha=0.35, linewidth=1.1)
ax1.yaxis.label.set_color('black')

# ---- Right axis : spray angle in degrees ----
ax2 = ax1.twinx()
h6, = ax2.plot(A_vec, alpha_n, 'b:', linewidth=2.2, label=r'$\alpha_n$  spray half-angle [right axis, °]')
ax2.set_ylabel(r'$\alpha$  (degrees)', fontsize=13, fontweight='bold')
ax2.set_ylim([0, 80])
ax2.set_yticks(np.arange(0, 90, 10))
ax2.tick_params(axis='y', labelsize=11)
ax2.yaxis.label.set_color('black')

# ---- Title ----
ax1.set_title(
    'Fig. 32 — Effects of Geometric Characteristic Parameter A\n'
    'on Injector Design and Flow Parameters',
    fontsize=12, fontweight='bold'
)

# ---- Legend ----
handles = [h1, h2, h3, h4, h5, h6]
ax1.legend(handles, [h.get_label() for h in handles],
           loc='center right', fontsize=10,
           frameon=True, edgecolor=(0.5, 0.5, 0.5))

# ---- Reference annotations ----
ax1.text(12, 0.97, r'$\bar{U}_{un}$', fontsize=11, color='b',            fontweight='bold')
ax1.text(12, 0.84, r'$\bar{h}$',      fontsize=11, color='k',            fontweight='bold')
ax1.text(12, 0.37, r'$\bar{U}_{an}$', fontsize=11, color='r',            fontweight='bold')
ax1.text(12, 0.20, r'$\phi$',         fontsize=11, color='m',            fontweight='bold')
ax1.text(12, 0.06, r'$\mu$',          fontsize=11, color=(0, 0.6, 0),    fontweight='bold')
ax2.text(12, 73,   r'$\alpha_n$',     fontsize=11, color='b',            fontweight='bold')

plt.tight_layout()
plt.savefig('swirl_injector_fig32.png', dpi=150, bbox_inches='tight')
print('\nFigure saved to: swirl_injector_fig32.png')


# ── 5.  Print reference values at A = 2, 4, 6, 8, 10 ──────────────────────
print(f"\n{'A':<5}  {'phi':<6}  {'mu':<6}  {'U_un':<6}  {'U_an':<6}  {'h_bar':<6}  {'alpha_n°':<8}")
print('-' * 55)
for A_ref in [2, 4, 6, 8, 10]:
    idx = np.argmin(np.abs(A_vec - A_ref))
    print(f"{A_ref:<5.0f}  {phi[idx]:<6.3f}  {mu[idx]:<6.3f}  "
          f"{U_un[idx]:<6.3f}  {U_an[idx]:<6.3f}  "
          f"{h_bar[idx]:<6.3f}  {alpha_n[idx]:<8.2f}")

plt.show()