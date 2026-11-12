# Pintle injector design
# Nitrous flows down the inside of the pintle and sprays outward (radially) through holes near its tip.
# Ethanol flows along the outside of the pintle through a ring of small axial holes (axial streams).
# The radial nitrous jets cross the axial ethanol streams, which is where they mix.
# (Swapped from ethanol-radial: the ethanol flow is too small for radial holes with good coverage.)
# Run cea_analysis.py, massIteration.py and heat_transfer.py first (or run_all.py).
import numpy as np
from CoolProp.CoolProp import PropsSI
import design as d
import two_phase as tp

psiToPa = 6894.76

# Inputs (from design_values.json)
Pc_Pa = d.get("Pc_psi") * psiToPa # chamber pressure (Pa)
mdot_ox = d.get("mdot_ox_kgps") # nitrous flow (kg/s)
mdot_fuel = d.get("mdot_fuel_kgps") # total ethanol flow (kg/s)
filmFraction = d.get("film_fraction") # fraction of ethanol going to the wall film
mdot_eth_core = mdot_fuel * (1 - filmFraction) # ethanol through the pintle holes (kg/s)
mdot_film = mdot_fuel * filmFraction # ethanol through the film holes (kg/s)
R_chamber = d.get("R_chamber_mm") / 1e3 # chamber radius (m)

# Injector pressure drops are a design choice. The pressure available at the injector minus the chosen
# injector drop is taken by the throttle valve, which gives it room to control the flow.
Cd = d.get("injector_Cd") # discharge coefficient

# Nitrous side: tank -> lines/valves -> throttle valve -> injector
P_ox_inj_in = (d.get("P_tank_psi") - d.get("feed_losses_psi")) * psiToPa # nitrous pressure available at injector (Pa)
dP_ox_available = P_ox_inj_in - Pc_Pa # most the nitrous injector could drop (Pa)
dP_ox = d.get("dP_ox_injector_psi") * psiToPa # chosen nitrous injector pressure drop (Pa)
rho_ox_liquid = d.get("rho_N2O") # liquid nitrous density in the tank (kg/m^3)
# The throttle valve drops the nitrous below its boiling pressure, so part of it boils on the way in:
# the injector gets a liquid + vapor mix (same enthalpy as the tank liquid). Two-phase (Dyer) flow, see two_phase.py.
h_ox = PropsSI('H', 'T', d.get("T_N2O_max_K"), 'P', P_ox_inj_in, 'NitrousOxide') # (J/kg)
P_ox_manifold = Pc_Pa + dP_ox # nitrous pressure inside the pintle (Pa)
rho_ox, x_ox = tp.state('NitrousOxide', P_ox_manifold, h_ox) # mix density and vapor fraction there

# Ethanol side: tank -> lines/valves -> cooling jacket -> throttle valve -> injector
P_eth_inj_in = d.get("P_ethanol_injector_in_psi") * psiToPa # ethanol pressure available at injector (Pa)
dP_eth_available = P_eth_inj_in - Pc_Pa # most the ethanol injector could drop (Pa)
dP_eth = d.get("dP_eth_injector_psi") * psiToPa # chosen ethanol injector pressure drop (Pa)
T_eth = d.get("T_coolant_out_K") # ethanol arrives hot from the jacket (K)
rho_eth = PropsSI('D', 'T', T_eth, 'P', Pc_Pa + dP_eth, 'Ethanol') # ethanol density at injector (kg/m^3)

print("-- Injector Inputs --")
print(f"Nitrous:   {mdot_ox*1e3:.2f} g/s, {x_ox*100:.0f} % vapor by mass after the throttle valve, mix density {rho_ox:.0f} kg/m^3 (liquid {rho_ox_liquid:.0f})")
print(f"           drop {dP_ox/psiToPa:.0f} psi ({dP_ox/Pc_Pa*100:.0f} % of Pc) of {dP_ox_available/psiToPa:.0f} psi available -> throttle valve takes {(dP_ox_available-dP_ox)/psiToPa:.0f} psi")
print(f"Ethanol:   {mdot_eth_core*1e3:.2f} g/s core + {mdot_film*1e3:.2f} g/s film, density {rho_eth:.0f} kg/m^3 at {T_eth:.0f} K")
print(f"           drop {dP_eth/psiToPa:.0f} psi ({dP_eth/Pc_Pa*100:.0f} % of Pc) of {dP_eth_available/psiToPa:.0f} psi available -> throttle valve takes {(dP_eth_available-dP_eth)/psiToPa:.0f} psi")
print(f"Chamber:   {2*R_chamber*1e3:.1f} mm diameter")
if dP_ox > dP_ox_available or dP_eth > dP_eth_available:
    print("WARNING: a chosen injector drop is bigger than the pressure available")

# Step 1: pintle diameter (chamber diameter / pintle diameter ~ 3-6)
D_chamber = 2 * R_chamber # m
D_pintle = D_chamber / d.get("pintle_ratio") # pintle outer diameter (m)

# Step 2: injection velocities, v = Cd * sqrt(2 * dP / rho); nitrous from its enthalpy drop (it boils as it expands)
v_ox = tp.jet_velocity('NitrousOxide', P_ox_manifold, h_ox, Pc_Pa, Cd) # nitrous, radial (m/s)
v_eth = Cd * np.sqrt(2 * dP_eth / rho_eth) # ethanol, axial (m/s)

# Step 3: flow areas. Ethanol: A = mdot / (rho * v) (the real hole area, since v already includes Cd).
# Nitrous: A = mdot / (Cd * G), with G the two-phase mass flux (Dyer). Liquid-only (SPI) shown for comparison.
A_ox = mdot_ox / (Cd * tp.mass_flux('NitrousOxide', P_ox_manifold, h_ox, Pc_Pa)) # total area of the nitrous radial holes (m^2)
A_ox_spi = mdot_ox / (Cd * np.sqrt(2 * rho_ox_liquid * dP_ox)) # what the first pass (all liquid) gave (m^2)
A_eth = mdot_eth_core / (rho_eth * v_eth) # total area of the ethanol axial holes (m^2)

# Step 4: ethanol axial holes (a ring of holes around the pintle base instead of a continuous gap,
# since a continuous gap for this little flow would be far too thin to make)
n_eth = int(d.get("n_eth_axial_holes")) # number of ethanol axial holes
d_eth = np.sqrt(4 * A_eth / (np.pi * n_eth)) # ethanol hole diameter (m)

# Step 5: total momentum ratio TMR = (radial momentum) / (axial momentum)
TMR = (mdot_ox * v_ox) / (mdot_eth_core * v_eth)

print("\n-- Pintle (Steps 1-5) --")
print(f"Pintle Diameter:             {D_pintle*1e3:.2f} mm (chamber/pintle = {d.get('pintle_ratio'):.1f})")
print(f"Nitrous Velocity (radial):   {v_ox:.1f} m/s")
print(f"Ethanol Velocity (axial):    {v_eth:.1f} m/s")
print(f"Nitrous Hole Area (total):   {A_ox*1e6:.3f} mm^2 (two-phase; all-liquid would be {A_ox_spi*1e6:.3f}, {A_ox/A_ox_spi:.1f}x smaller)")
print(f"Ethanol Hole Area (total):   {A_eth*1e6:.3f} mm^2")
print(f"Ethanol Axial Holes:         {n_eth} x {d_eth*1e3:.3f} mm")
print(f"Total Momentum Ratio (TMR):  {TMR:.2f} (ideal ~1, workable ~0.5-2.5)")
if TMR > 2.5:
    print("WARNING: TMR too high: the nitrous jets punch through the ethanol streams toward the wall")

# Step 6: nitrous radial holes (number, diameter) and blockage factor
# Blockage factor = fraction of the pintle circumference covered by holes = (N * d) / (pi * D_pintle), aim ~0.3-0.7
n_ox = int(d.get("n_ox_radial_holes")) # number of nitrous radial holes
d_ox = np.sqrt(4 * A_ox / (np.pi * n_ox)) # nitrous hole diameter (m)
blockage = n_ox * d_ox / (np.pi * D_pintle) # blockage factor
ligament = np.pi * D_pintle / n_ox - d_ox # solid metal left between neighbouring holes (m)
r_bore = d.get("pintle_bore_d_m") / 2 # pintle nitrous bore radius (m)
ligament_in = 2 * r_bore * np.sin(np.pi / n_ox) - d_ox # same, where the holes break into the bore (m)

print("\n-- Pintle (Step 6) --")
print(f"Nitrous Radial Holes:        {n_ox} x {d_ox*1e3:.3f} mm")
print(f"Blockage Factor:             {blockage:.2f} (aim 0.3-0.7)")
print(f"Metal Between Holes:         {ligament*1e3:.2f} mm outside, {ligament_in*1e3:.2f} mm at the bore")
if ligament_in < 0.3e-3:
    print("WARNING: neighbouring nitrous holes (nearly) run into each other inside the pintle")
if blockage < 0.3 or blockage > 0.7:
    print("WARNING: blockage factor outside 0.3-0.7")

d.save_many({
    "D_pintle_mm": (D_pintle*1e3, "mm"),
    "v_ox_injector_mps": (v_ox, "m/s"),
    "v_eth_injector_mps": (v_eth, "m/s"),
    "A_ox_holes_mm2": (A_ox*1e6, "mm^2"),
    "A_eth_holes_mm2": (A_eth*1e6, "mm^2"),
    "d_eth_axial_hole_mm": (d_eth*1e3, "mm"),
    "d_ox_radial_hole_mm": (d_ox*1e3, "mm"),
    "TMR": (TMR, "-"),
    "blockage_factor": (blockage, "-"),
    "rho_ox_manifold": (rho_ox, "kg/m^3"),
    "x_ox_manifold": (x_ox, "-"),
}, "injector.py")

# Step 7: film holes (ring at the outer edge of the face)
# The film ethanol comes from the same ethanol manifold as the core (manifold pressure = Pc + dP_eth).
# Path: manifold -> one restrictor hole -> small film manifold -> ring of film holes -> along the wall.
# The film holes only take a small drop so they can be bigger; the single restrictor takes the rest.
n_film = int(d.get("n_film_holes")) # number of film holes
dP_film_holes = d.get("dP_film_holes_psi") * psiToPa # drop across the film holes (Pa)
dP_restrictor = dP_eth - dP_film_holes # drop across the restrictor (Pa)
rho_film = rho_eth # same hot ethanol as the core (kg/m^3)

v_film = Cd * np.sqrt(2 * dP_film_holes / rho_film) # film jet velocity (m/s)
A_film = mdot_film / (rho_film * v_film) # total film hole area (m^2)
d_film = np.sqrt(4 * A_film / (np.pi * n_film)) # film hole diameter (m)
R_film = R_chamber - d.get("film_hole_wall_gap_m") # radius of the film hole ring (m)
film_spacing = 2 * np.pi * R_film / n_film # distance between neighbouring film holes along the ring (m)

v_restrictor = Cd * np.sqrt(2 * dP_restrictor / rho_film) # velocity through the restrictor (m/s)
A_restrictor = mdot_film / (rho_film * v_restrictor) # restrictor area (m^2)
d_restrictor = np.sqrt(4 * A_restrictor / np.pi) # single restrictor hole diameter (m)

print("\n-- Film Holes (Step 7) --")
print(f"Film Flow:                   {mdot_film*1e3:.2f} g/s ({filmFraction*100:.0f} % of ethanol)")
print(f"Film Holes:                  {n_film} x {d_film*1e3:.3f} mm, ring radius {R_film*1e3:.2f} mm, {film_spacing*1e3:.1f} mm apart")
print(f"Film Jet Velocity:           {v_film:.1f} m/s (drop {dP_film_holes/psiToPa:.0f} psi)")
print(f"Restrictor:                  1 x {d_restrictor*1e3:.3f} mm (drop {dP_restrictor/psiToPa:.0f} psi)")
if dP_restrictor <= 0:
    print("WARNING: film hole drop is bigger than the whole ethanol injector drop")

d.save_many({
    "d_film_hole_mm": (d_film*1e3, "mm"),
    "R_film_ring_mm": (R_film*1e3, "mm"),
    "film_hole_spacing_mm": (film_spacing*1e3, "mm"),
    "d_film_restrictor_mm": (d_restrictor*1e3, "mm"),
}, "injector.py")

# Step 8: skip distance (how far the pintle sticks out) ~ 1 x D_pintle
# The axial ethanol streams travel this far along the pintle before the radial nitrous jets cross them.
# The nitrous holes sit at the skip distance; the pintle's bore stops at the holes and a closed end wall
# (thickness checked in injector_hardware.py) caps the tip.
L_skip = d.get("skip_ratio") * D_pintle # face to the centre of the nitrous holes (m)
L_tip = L_skip + d_ox / 2 + d.get("pintle_tip_wall_m") # face to the pintle tip (m)
radial_travel = R_chamber - D_pintle / 2 # how far the nitrous jets travel from the pintle to the wall (m)
L_cyl = d.get("L_cyl_mm") / 1e3 # straight chamber length (m)

print("\n-- Skip Distance (Step 8) --")
print(f"Skip Distance:               {L_skip*1e3:.2f} mm (face to nitrous holes, {d.get('skip_ratio'):.1f} x pintle diameter)")
print(f"Pintle Tip From Face:        {L_tip*1e3:.2f} mm")
print(f"Nitrous Jet Travel to Wall:  {radial_travel*1e3:.2f} mm")
print(f"Straight Chamber Length:     {L_cyl*1e3:.2f} mm (pintle tip uses {L_tip/L_cyl*100:.0f} % of it)")
if L_tip > 0.5 * L_cyl:
    print("WARNING: pintle reaches past half the straight chamber, less room to burn before the throat")

d.save_many({
    "L_skip_mm": (L_skip*1e3, "mm"),
    "L_pintle_tip_mm": (L_tip*1e3, "mm"),
}, "injector.py")
