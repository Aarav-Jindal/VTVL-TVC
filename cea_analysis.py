import numpy as np
import cea
import design as d

# Inputs (from design_values.json)
Pc = d.get("Pc_psi") #chamber pressure (psia)
Pe = d.get("Pe_psi") #exit pressure (psia)
Pi_p = Pc/Pe #pressure ratio
Pc = cea.units.psi_to_bar(Pc) #convert to bar
Pe = cea.units.psi_to_bar(Pe) #convert to bar
mixtureRatios = np.arange(2,8,0.25) #mixture ratios to evaluate

#Reactants
reac_names = ["N2O", "C2H5OH(L)"] #N2O and Ethanol
fuel_weights = np.array([0.0, 1.0]) # Weight fractions of fuel
oxidant_weights = np.array([1.0, 0.0]) # Weight fractions of oxidant
T_reactant = np.array([298.15, 298.15])  # Reactant temperatures (K)

#CEA Objects
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)
solver = cea.RocketSolver(prod, reactants=reac) # performance (Isp, c*, Ae/At)
soln = cea.RocketSolution(solver)
# Separate solver for transport properties only: with transport=True this CEA version
# returns NaN/wrong nozzle results, but the chamber values (index 0) are fine
solverTransport = cea.RocketSolver(prod, reactants=reac, transport=True)
solnTransport = cea.RocketSolution(solverTransport)

#Empty lists to store results
Isp = []
c_star = []
Tc = []
gamma = []
Ae_At = []
viscosity = []
cp_fr = []
Pr_fr = []

# Loop through mixture ratios and solve
for of_ratio in mixtureRatios:
    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R
    solver.solve(soln, weights, Pc, Pi_p, hc=hc, n_frz=2)
    Isp.append(soln.Isp[-1]/9.80665) #convert to seconds
    c_star.append(soln.c_star[-1]*3.28084) #convert to ft/s
    Tc.append(soln.T[0])
    gamma.append(soln.gamma_s[0])
    Ae_At.append(soln.ae_at[-1])
    solverTransport.solve(solnTransport, weights, Pc, Pi_p, hc=hc, n_frz=2)
    viscosity.append(solnTransport.viscosity[0])
    cp_fr.append(solnTransport.cp_fr[0])
    Pr_fr.append(solnTransport.Pr_fr[0])


# Print results
print("Mixture Ratio\tIsp (s)\tc* (ft/s)\tTc (K)\tgamma\tAe/At")
for i in range(len(mixtureRatios)):
    print(f"{mixtureRatios[i]:.2f}\t\t{Isp[i]:.2f}\t\t{c_star[i]:.2f}\t\t{Tc[i]:.2f}\t{gamma[i]:.4f}\t{Ae_At[i]:.4f}")

# Print best mixture ratio
best_index = np.argmax(Isp)
print(f"\nBest mixture ratio for maximum Isp: {mixtureRatios[best_index]:.2f} with Isp = {Isp[best_index]:.2f} s")

# Design Parameters (from design_values.json)
F = d.get("thrust_N") # thrust (N), fixed input in design_values.json
mixRatio_chosen = d.get("mixture_ratio")
effeciency = d.get("efficiency")
g0 = 9.80665 # m/s^2

# 1. Find the row for your mixture ratio
i = np.argmin(np.abs(mixtureRatios - mixRatio_chosen))
Isp_chosen = Isp[i]
cstar_chosen = c_star[i]/3.28084 # convert back to m/s
mu_chosen = viscosity[i] * 1e-4   # millipoise -> Pa·s
cp_chosen = cp_fr[i] * 1000       # kJ/(kg·K) -> J/(kg·K)
Pr_chosen = Pr_fr[i]
gamma_chosen = gamma[i]

## applying effeciency factor to Isp and cstar
Isp_real = Isp_chosen * effeciency
cstar_real = cstar_chosen * effeciency

Pc_Pa = Pc * 1e5 # convert bar to Pa

# 2. Calculate mass flow rate
mdot = F / (Isp_real * g0) # kg/s

# Split into oxidizer and fuel mass flow rates
mdot_ox = mdot * (mixRatio_chosen / (1 + mixRatio_chosen)) # kg/s
mdot_fuel = mdot * (1 / (1 + mixRatio_chosen))

# Throat area and diameter
At = (mdot * cstar_real) / Pc_Pa # m^2
Dt = np.sqrt((4 * At) / np.pi) # m

# Exit area and diameter
Ae = At * Ae_At[i] # m^2
De = np.sqrt((4 * Ae) / np.pi) # m

# Chamber Geometry
L_star = d.get("L_star_m") # characteristic length (m)
contraction_ratio = d.get("contraction_ratio") # contraction ratio
V_chamber = L_star * At # m^3
A_chamber = At * contraction_ratio # m^2
L_chamber = V_chamber / A_chamber # m

# Nozzle Geometry
divergence_angle = np.arange(5, 30, 1) # divergence angles to evaluate (degrees)
lossFactor = (1+np.cos(np.radians(divergence_angle)))/2 # lambda factor for each angle
L = (De - Dt) / (2 * np.tan(np.radians(divergence_angle))) # nozzle length for each angle

# Plot Lambda vs Divergence Angle
# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 6))
# plt.plot(divergence_angle, lossFactor, 'o-')
# plt.xlabel('Divergence Angle (degrees)')
# plt.ylabel('Loss Factor')
# plt.title('Loss Factor vs Divergence Angle')
# plt.grid(True)
# plt.show()

# Nozzle Geometry for chosen divergence angle
chosen_angle = d.get("diverging_angle_deg") # degrees
convering_angle = d.get("converging_angle_deg") # degrees
R_chamber = np.sqrt(A_chamber / np.pi) # chamber radius (m)
R_throat = Dt / 2 # throat radius (m)
R_exit = De / 2 # exit radius (m)
L_nozzle = (De - Dt) / (2 * np.tan(np.radians(chosen_angle))) # nozzle length (m)
L_converging = (R_chamber-R_throat) / np.tan(np.radians(convering_angle)) # converging section length (m)
V_converging = (1/3) * np.pi * L_converging * (R_chamber**2 + R_chamber*R_throat + R_throat**2) # converging section volume (m^3)
L_cyl = (V_chamber - V_converging) / A_chamber # cylindrical section length (m)

# Print design parameters
L_total = L_cyl + L_converging + L_nozzle # injector face to nozzle exit (m)
lambda_chosen = (1 + np.cos(np.radians(chosen_angle))) / 2 # divergence efficiency for chosen angle

print(f"\nDesign Parameters for Mixture Ratio {mixRatio_chosen}:")
print("\n-- Performance --")
print(f"Isp (real):                  {Isp_real:.2f} s")
print(f"c* (real):                   {cstar_real:.1f} m/s")
print(f"Mass Flow Rate:              {mdot*1e3:.2f} g/s")
print(f"Oxidizer Mass Flow Rate:     {mdot_ox*1e3:.2f} g/s")
print(f"Fuel Mass Flow Rate:         {mdot_fuel*1e3:.2f} g/s")

print("\n-- Areas and volume --")
print(f"Throat Area:                 {At*1e6:.2f} mm^2")
print(f"Exit Area:                   {Ae*1e6:.2f} mm^2")
print(f"Chamber Area:                {A_chamber*1e6:.2f} mm^2")
print(f"Chamber Volume:              {V_chamber*1e6:.2f} cm^3")
print(f"Converging Section Volume:   {V_converging*1e6:.2f} cm^3")

print("\n-- Radii (CAD) --")
print(f"Chamber Radius:              {R_chamber*1e3:.4f} mm")
print(f"Throat Radius:               {R_throat*1e3:.4f} mm")
print(f"Exit Radius:                 {R_exit*1e3:.4f} mm")

print("\n-- Lengths --")
print(f"Cylindrical Section Length:  {L_cyl*1e3:.4f} mm")
print(f"Converging Section Length:   {L_converging*1e3:.2f} mm")
print(f"Nozzle Length:               {L_nozzle*1e3:.2f} mm")
print(f"Total Length (face to exit): {L_total*1e3:.2f} mm")

print("\n-- Angles --")
print(f"Converging Angle:            {convering_angle} deg")
print(f"Divergence Angle:            {chosen_angle} deg")
print(f"Divergence Efficiency:       {lambda_chosen:.4f}")

print("\n-- Gas Properties (Chamber) --")
print(f"Viscosity:                   {mu_chosen:.3e} Pa·s")
print(f"Specific Heat (cp):          {cp_chosen:.2f} J/(kg·K)")
print(f"Prandtl Number:              {Pr_chosen:.4f}")
print(f"Gamma:                       {gamma_chosen:.4f}")

# Save results for the other scripts
d.save_many({
    "Isp_ideal_s": (Isp_chosen, "s"),
    "Isp_real_s": (Isp_real, "s"),
    "cstar_real_mps": (cstar_real, "m/s"),
    "Tc_K": (Tc[i], "K"),
    "gamma": (gamma_chosen, "-"),
    "mu_gas_Pa_s": (mu_chosen, "Pa*s"),
    "cp_gas_J_kgK": (cp_chosen, "J/(kg*K)"),
    "Pr_gas": (Pr_chosen, "-"),
    "Ae_At": (Ae_At[i], "-"),
    "mdot_total_kgps": (mdot, "kg/s"),
    "mdot_ox_kgps": (mdot_ox, "kg/s"),
    "mdot_fuel_kgps": (mdot_fuel, "kg/s"),
    "A_throat_mm2": (At*1e6, "mm^2"),
    "R_chamber_mm": (R_chamber*1e3, "mm"),
    "R_throat_mm": (R_throat*1e3, "mm"),
    "R_exit_mm": (R_exit*1e3, "mm"),
    "L_cyl_mm": (L_cyl*1e3, "mm"),
    "L_converging_mm": (L_converging*1e3, "mm"),
    "L_nozzle_mm": (L_nozzle*1e3, "mm"),
    "L_total_mm": (L_total*1e3, "mm"),
    "V_chamber_cm3": (V_chamber*1e6, "cm^3"),
}, "cea_analysis.py")

#Plot results
import matplotlib.pyplot as plt

## Figure 1: Isp vs Mixture Ratio
# plt.figure(figsize=(10, 6))
# plt.plot(mixtureRatios, Isp, 'o-')
# plt.xlabel('Mixture Ratio')
# plt.ylabel('Isp (s)')
# plt.title('Performance of Different Mixture Ratios')
# plt.grid(True)
# plt.show()

# ## Figure 2: Tc vs Mixture Ratio
# plt.figure(figsize=(10, 6))
# plt.plot(mixtureRatios, Tc, 'o-')
# plt.xlabel('Mixture Ratio')
# plt.ylabel('Tc (K)')
# plt.title('Combustion Temperature vs Mixture Ratio')
# plt.grid(True)
# plt.show()
