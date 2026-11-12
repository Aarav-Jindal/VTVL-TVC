import csv, os
import numpy as np
import design as d

g = 9.80665 #m/s^2

# Hardware before margin: the parts list (parts_list.csv, known vs estimated) + the engine parts from the CAD
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "parts_list.csv"), newline="", encoding="utf-8") as f:
    parts = list(csv.DictReader(line for line in f if not line.startswith("#")))
listKnown = sum(int(p["qty"]) * float(p["each_g"]) for p in parts if p["group"] == "known") / 1000 #kg
hardwareEst = sum(int(p["qty"]) * float(p["each_g"]) for p in parts if p["group"] == "est") / 1000 #kg, estimates and placeholders
hardwareKnown = listKnown + d.get("m_engine_cad_kg") + d.get("m_gimbal_cad_kg") #kg, real masses: catalog / calculated parts + the engine and gimbal CAD

# Inputs (from design_values.json)
marginKnown = d.get("margin_known") # on known hardware and the tanks
marginEst = d.get("margin_est") # on estimated hardware and the helium bottle
mixtureRatio = d.get("mixture_ratio")
rho_N2O = d.get("rho_N2O")
rho_ethanol = d.get("rho_ethanol")
ullageFactor = d.get("ullage_factor")
P_tank = d.get("P_tank_psi") #psi
SF = d.get("tank_safety_factor")
sigma = d.get("tank_yield_Pa") #Pa
rho_Al = 2700 #kg/m^3
bulkheadMass = d.get("bulkhead_mass_kg") #kg
bulkheadLengthAllowance = d.get("bulkhead_length_m") #m
R_He = 2077 #J/kgK
T = 298 #K
He_margin = d.get("He_margin")
P_full = d.get("He_bottle_full_psi") #psi
P_min = P_tank + d.get("regulator_dropout_psi") #psi, regulator needs this much above its output to keep regulating

r_out = d.get("tank_outer_radius_m") #m, 4" OD tube
wall = d.get("tank_wall_m") #m, 0.125" wall
r_in = r_out - wall #m

# Body tube: its length is the sections that don't change + both tanks (same stack as vehicle_layout.py), so it grows with them
bodyKgPerM = d.get("body_tube_kg_per_m") #kg/m
L_fixed = (d.get("engine_bay_len_m") + d.get("He_valve_stack_len_m") + d.get("He_bottle_len_m") + d.get("avionics_bay_len_m")
           + 3 * d.get("section_gap_m") + d.get("n2o_eth_gap_m") + 2 * bulkheadLengthAllowance) #m, height of the whole stack without the tanks
# The tube stops just below the gimbal mount: the engine hangs out the bottom, so its outlet fitting has no wall to hit.
# The extra length below the mount gives the mount's screw holes enough fiberglass between them and the tube's end
x_bodyBottom = d.get("L_gas_total_mm") / 1000 + d.get("gimbal_pivot_above_face_m") - d.get("gimbal_mount_h_m") / 2 - d.get("body_tube_below_mount_m") #m above the nozzle exit

# Propellant fraction from the rocket equation, using the latest Isp from cea_analysis.py
deltaV = d.get("delta_v_mps") # m/s, from the flight profile
Isp_real = d.get("Isp_real_s") # s
massRatio = np.exp(deltaV / (Isp_real * g))
propellantFraction = 1 - 1 / massRatio

bottles = [ # (name, size in ci, mass in kg), smallest first
    ("HK Army HSTL 98/4500", 98, 0.93),
    ("114/4500 (mass est.)", 114, 1.10),
    ("First Strike Hero 2.0 124/4500 (mass est.)", 124, 1.25), # no published mass; scaled from the 98 ci bottle, rounded up
    ("3L 183/4500 carbon fiber (VXB listing, CE not DOT)", 183, 2.0), # seller-listed mass, not a spec sheet
    ("Luxfer L45R 285/4500 (DOT, 5.4 x 18.1 in)", 285, 2.68), # Luxfer spec sheet 5.9 lb, no valve or adapter
    ("Luxfer ECLIPSE PRO L66P 417/4500 (6.1 x 20.3 in)", 417, 3.40), # Luxfer spec sheet 7.5 lb
    ("Luxfer ECLIPSE PRO L88D 548/4500 (6.9 x 21.5 in)", 548, 4.58), # Luxfer spec sheet 10.1 lb
    ("2 x Luxfer ECLIPSE PRO L66P 834/4500 in line", 834, 6.80),
]

thrust = d.get("thrust_N") #N, fixed engine thrust (input, not set by this loop)
minThrustToWeight = d.get("thrust_to_weight") # lowest T/W allowed at liftoff before the engine needs resizing
throttleMin = d.get("throttle_min") # lowest throttle the engine can run at

m_sizingFrozen = d.get("m_sizing_frozen_kg") #kg, the mass the propellant and tanks are sized for
m_wet = 15.0 #kg, starting guess
maxIterations = 50 # set to 1 to compare with a hand calculation

# Unit conversions
psiToPa = 6894.76
ciToM3 = 1.6387e-5
P_tank_Pa = P_tank * psiToPa
P_full_Pa = P_full * psiToPa
P_min_Pa = P_min * psiToPa

# Check the chosen tube wall is thick enough (hoop stress)
t_required = P_tank_Pa * SF * r_out / sigma
if wall < P_tank_Pa * d.get("tank_wall_SF_accepted") * r_out / sigma: # SF 1.99 on this wall accepted 2026-11-05
    print(f"WARNING: wall {wall*1e3:.2f} mm is thinner than required {t_required*1e3:.2f} mm")

A_in = np.pi * r_in**2 # tank inner cross-section (m^2)
A_metal = np.pi * (r_out**2 - r_in**2) # tube wall cross-section (m^2)

converged = False
for iteration in range(1, maxIterations + 1):
    # Propellant. The tanks are frozen at this mass (decided 2026-11-06): a new mass in the parts list changes T/W and apogee, not the tanks, so the CAD stays in step
    m_prop = propellantFraction * m_sizingFrozen
    m_N2O = m_prop * mixtureRatio / (1 + mixtureRatio)
    m_ethanol = m_prop - m_N2O

    # Tank volumes and lengths
    V_N2O = m_N2O / rho_N2O * ullageFactor # m^3
    V_ethanol = m_ethanol / rho_ethanol * ullageFactor # m^3
    L_N2O = V_N2O / A_in # m
    L_ethanol = V_ethanol / A_in # m

    # Tank masses (tube + 2 bulkheads each)
    m_tank_N2O = A_metal * (L_N2O + bulkheadLengthAllowance) * rho_Al + 2 * bulkheadMass
    m_tank_ethanol = A_metal * (L_ethanol + bulkheadLengthAllowance) * rho_Al + 2 * bulkheadMass
    L_body = L_fixed + L_N2O + L_ethanol - x_bodyBottom #m
    m_body = bodyKgPerM * L_body #kg

    # Helium needed to fill both tanks at tank pressure
    m_He = P_tank_Pa * (V_N2O + V_ethanol) / (R_He * T) * He_margin

    # Bottle size from the usable pressure range, then pick the smallest bottle that fits
    V_bottle = m_He * R_He * T / (P_full_Pa - P_min_Pa) # m^3
    bottleNeeded = V_bottle / ciToM3 # ci
    bottle = next((b for b in bottles if b[1] >= bottleNeeded), None)
    if bottle is None:
        print(f"WARNING: need a {bottleNeeded:.1f} ci bottle, none in the list is big enough")
        break
    bottleName, bottleSize, bottleMass = bottle

    # Total
    m_dry = (hardwareKnown + m_tank_N2O + m_tank_ethanol) * marginKnown + (hardwareEst + m_body + bottleMass) * marginEst
    m_wet_new = m_dry + m_prop + m_He

    converged = abs(m_wet_new - m_wet) < 0.001
    m_wet = m_wet_new
    if converged:
        break

# Check the fixed engine still fits this vehicle
thrustToWeight = thrust / (m_wet * g) # at liftoff, full tanks
m_land_light = hardwareKnown + hardwareEst + m_body + m_tank_N2O + m_tank_ethanol + bottleMass + m_He # empty vehicle with none of the margin used (lightest landing)
throttleLanding = m_land_light * g / thrust # throttle needed to hover it just before touchdown

# Print results
print(f"\nPropellant fraction:     {propellantFraction:.4f} (delta-v {deltaV} m/s, Isp {Isp_real:.2f} s)")
print(f"\nIterations: {iteration} ({'converged' if converged else 'NOT converged'})")
print("\n-- Masses --")
print(f"Tanks frozen at:         {m_sizingFrozen:.3f} kg. Real delta-v at this wet mass: {Isp_real*g*np.log(m_wet/(m_wet - m_prop)):.1f} m/s (the sizing assumed {deltaV})")
print(f"Wet mass:                {m_wet:.3f} kg")
print(f"Dry mass (with margin):  {m_dry:.3f} kg")
print(f"Propellant:              {m_prop:.3f} kg")
print(f"  N2O:                   {m_N2O:.3f} kg")
print(f"  Ethanol:               {m_ethanol:.3f} kg")
print(f"Helium:                  {m_He*1e3:.1f} g")
print(f"Hardware (no margin):    {listKnown:.3f} kg known + {d.get('m_engine_cad_kg'):.3f} kg engine CAD + {d.get('m_gimbal_cad_kg'):.3f} kg gimbal CAD + {hardwareEst:.3f} kg estimated (parts_list.csv)")
print(f"Body tube:               {L_body:.2f} m, {m_body:.2f} kg")

print("\n-- Tanks --")
print(f"Wall required / chosen:  {t_required*1e3:.2f} / {wall*1e3:.2f} mm")
print(f"N2O tank:                {L_N2O*100:.1f} cm long, {m_tank_N2O:.3f} kg")
print(f"Ethanol tank:            {L_ethanol*100:.1f} cm long, {m_tank_ethanol:.3f} kg")

print("\n-- Helium bottle --")
print(f"Needed:                  {bottleNeeded:.1f} ci")
print(f"Chosen:                  {bottleName}, {bottleMass:.2f} kg")

print("\n-- Engine fit (thrust is a fixed input) --")
print(f"Thrust:                  {thrust:.1f} N = {thrust/g:.2f} kgf")
print(f"T/W at liftoff:          {thrustToWeight:.2f} (min {minThrustToWeight})")
print(f"Landing throttle:        {throttleLanding*100:.0f} % with none of the margin used (min {throttleMin*100:.0f} %)")
if thrustToWeight < minThrustToWeight:
    print("WARNING: T/W at liftoff is below the minimum, the engine needs more thrust")
if throttleLanding < throttleMin:
    print("WARNING: the engine can't throttle low enough to land the lightest vehicle")

# Save results for the other scripts (only a converged result)
if converged:
    d.save_many({
        "propellant_fraction": (propellantFraction, "-"),
        "m_wet_kg": (m_wet, "kg"),
        "m_dry_kg": (m_dry, "kg"),
        "m_prop_kg": (m_prop, "kg"),
        "m_N2O_kg": (m_N2O, "kg"),
        "m_ethanol_kg": (m_ethanol, "kg"),
        "m_He_kg": (m_He, "kg"),
        "L_tank_N2O_m": (L_N2O, "m"),
        "L_tank_ethanol_m": (L_ethanol, "m"),
        "m_tank_N2O_kg": (m_tank_N2O, "kg"),
        "m_tank_ethanol_kg": (m_tank_ethanol, "kg"),
        "He_bottle_needed_ci": (bottleNeeded, "ci"),
        "He_bottle": (bottleName, "-"),
        "He_bottle_mass_kg": (bottleMass, "kg"),
        "thrust_to_weight_liftoff": (thrustToWeight, "-"),
        "throttle_landing_light": (throttleLanding, "-"),
        "hardware_known_kg": (listKnown, "kg"),
        "hardware_est_kg": (hardwareEst, "kg"),
        "L_body_tube_m": (L_body, "m"),
        "body_tube_bottom_m": (x_bodyBottom, "m"),
        "m_body_tube_kg": (m_body, "kg"),
    }, "massIteration.py")
else:
    print("\nNOT saved: the mass loop did not converge")
