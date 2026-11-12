"""
Flow of a liquid that can boil (flash) as it passes through a hole or valve.
Nitrous always does this in the injector; hot ethanol does it at low throttle.

Dyer model (Dyer et al. 2007): a blend of two limits, weighted by how fast bubbles can grow:
  SPI  single-phase incompressible: the liquid doesn't boil inside the hole (too much flow)
  HEM  homogeneous equilibrium: the liquid boils instantly as the pressure falls (too little flow)
If the fluid is already part vapor when it reaches the hole (it boiled in a valve upstream),
bubbles are already there, so HEM alone is used.

    import two_phase as tp
    G = tp.mass_flux("NitrousOxide", P1, h1, P2)   # kg/(m^2 s); mdot = Cd * area * G
"""
import numpy as np
from CoolProp.CoolProp import PropsSI


def state(fluid, P, h):
    """Density (kg/m^3) and vapor mass fraction (0 = all liquid) at pressure P (Pa) and enthalpy h (J/kg)."""
    q = PropsSI('Q', 'P', P, 'H', h, fluid)
    return PropsSI('D', 'P', P, 'H', h, fluid), (q if 0 < q <= 1 else 0.0)


def vapor_pressure(fluid, P, h):
    """Pressure the liquid starts boiling at, for the liquid at P (Pa) and h (J/kg)."""
    T = PropsSI('T', 'P', P, 'H', h, fluid)
    return PropsSI('P', 'T', T, 'Q', 0, fluid) if T < PropsSI('Tcrit', fluid) else np.inf


def hem_flux(fluid, P1, h1, P2):
    """HEM mass flux (kg/(m^2 s)): expands with no losses from P1 and boils instantly.
    If the flux peaks above P2 the hole is choked (lowering P2 further adds no flow)."""
    s1 = PropsSI('S', 'P', P1, 'H', h1, fluid)
    G = 0.0
    for P in np.linspace(P2, P1, 40, endpoint=False):
        h = PropsSI('H', 'P', P, 'S', s1, fluid)
        G = max(G, PropsSI('D', 'P', P, 'S', s1, fluid) * np.sqrt(2 * max(h1 - h, 0)))
    return G


def mass_flux(fluid, P1, h1, P2):
    """Mass flux through a hole (kg/(m^2 s)) from P1 (Pa, enthalpy h1 J/kg) into P2 (Pa). mdot = Cd * area * flux."""
    rho1, x1 = state(fluid, P1, h1)
    if x1 > 0: # already part vapor: no delay to form bubbles
        return hem_flux(fluid, P1, h1, P2)
    G_spi = np.sqrt(2 * rho1 * (P1 - P2))
    Pv = min(vapor_pressure(fluid, P1, h1), P1)
    if P2 >= Pv: # never drops low enough to boil
        return G_spi
    k = np.sqrt((P1 - P2) / (Pv - P2)) # time for bubbles to grow / time spent in the hole
    return k / (1 + k) * G_spi + 1 / (1 + k) * hem_flux(fluid, P1, h1, P2)


def jet_velocity(fluid, P1, h1, P2, Cd):
    """Jet velocity once fully expanded to P2 (m/s), from the enthalpy drop. For a liquid this is Cd * sqrt(2 dP / rho)."""
    h2 = PropsSI('H', 'P', P2, 'S', PropsSI('S', 'P', P1, 'H', h1, fluid), fluid)
    return Cd * np.sqrt(2 * (h1 - h2))
