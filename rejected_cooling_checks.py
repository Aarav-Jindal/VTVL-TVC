# Cooling options that were checked and rejected, and why.
# 1. Annular-gap jacket at the throat: too little ethanol flow, so the coolant is too slow (h_c far too low).
# 2. Copper heat-sink throat: the hot surface passes the copper limit within seconds, far short of a ~31 s burn.
# Chosen instead: lower chamber pressure (300 psi) + milled regen channels + 20% film (see heat_transfer.py).
# Run heat_transfer.py first: this uses its throat heat flux.
import numpy as np
from CoolProp.CoolProp import PropsSI
import design as d

psiToPa = 6894.76
R_throat = d.get("R_throat_mm") / 1e3 # m
mdot_fuel = d.get("mdot_fuel_kgps") # kg/s
T_wg = d.get("T_wg_K") # K
q_throat_film = d.get("q_throat_film_W_m2") # film-cooled heat flux at the throat (W/m^2)
h_g_throat = d.get("h_g_throat_W_m2K") # Bartz h_g at the throat (W/(m^2*K))

# Cooling jacket feasibility check at the throat (annular gap between copper liner and outer sleeve)
t_liner = 1e-3 # copper liner thickness (m), placeholder
k_copper = 390 # copper thermal conductivity (W/(m*K))
T_coolant = 350 # ethanol temperature in the jacket at the throat (K), placeholder
P_jacket = (d.get("P_tank_psi") - d.get("feed_losses_psi")) * psiToPa # ethanol pressure entering the jacket (Pa)
gap = 0.3e-3 # annular coolant gap (m), placeholder

# Step 1: coolant-side heat transfer coefficient the throat needs
dT_liner = q_throat_film * t_liner / k_copper # temperature drop through the liner (K)
T_wc = T_wg - dT_liner # coolant-side wall temperature (K)
h_c_needed = q_throat_film / (T_wc - T_coolant) # W/(m^2*K)

# Step 2: coolant-side heat transfer coefficient the ethanol flow actually gives
rho_cool = PropsSI('D', 'T', T_coolant, 'P', P_jacket, 'Ethanol') # density (kg/m^3)
mu_cool = PropsSI('V', 'T', T_coolant, 'P', P_jacket, 'Ethanol') # viscosity (Pa*s)
k_cool = PropsSI('L', 'T', T_coolant, 'P', P_jacket, 'Ethanol') # thermal conductivity (W/(m*K))
Pr_cool = PropsSI('PRANDTL', 'T', T_coolant, 'P', P_jacket, 'Ethanol') # Prandtl number

A_flow = 2 * np.pi * (R_throat + t_liner) * gap # flow area of the annular gap (m^2)
v_cool = mdot_fuel / (rho_cool * A_flow) # coolant velocity (m/s), all the ethanol goes through the jacket
D_h = 2 * gap # hydraulic diameter of a thin annular gap (m)
Re_cool = rho_cool * v_cool * D_h / mu_cool # Reynolds number
Nu_cool = 0.023 * Re_cool**0.8 * Pr_cool**0.4 # Nusselt number (Dittus-Boelter)
h_c_achieved = Nu_cool * k_cool / D_h # W/(m^2*K)

print("\n-- Cooling Jacket Feasibility (Throat) --")
print(f"Throat Heat Flux (film):     {q_throat_film/1e6:.2f} MW/m^2")
print(f"Liner Temperature Drop:      {dT_liner:.1f} K  ->  coolant-side wall {T_wc:.1f} K")
print(f"Coolant Velocity:            {v_cool:.2f} m/s  (gap {gap*1e3:.2f} mm, Re {Re_cool:.0f})")
print(f"h_c Needed:                  {h_c_needed/1e3:.1f} kW/(m^2*K)")
print(f"h_c Achieved:                {h_c_achieved/1e3:.1f} kW/(m^2*K)")
print(f"Feasible?                    {'YES' if h_c_achieved >= h_c_needed else 'NO'}")

# Copper heat-sink throat: transient heat soak (1D radial conduction through a copper ring)
# Inner surface heated by the hot gas (Bartz + film), outer surface cooled by the regen ethanol.
# Axial conduction is ignored (conservative: in reality heat also spreads along the engine).
r_block_out = 15e-3 # outer radius of the copper throat block (m), placeholder
burn_time = 12.88 + 12.88 + 5 # ascent + landing + hover burn, run continuously at full thrust (s), conservative
T_copper_limit = 800 # max hot-surface temperature for copper (K)
rho_cu = 8940 # copper density (kg/m^3)
cp_cu = 385 # copper specific heat (J/(kg*K))
alpha_cu = k_copper / (rho_cu * cp_cu) # thermal diffusivity (m^2/s)

# Hot side: Bartz h_g at the throat and the film temperature the wall sees there
T_seen_throat = q_throat_film / h_g_throat + T_wg

# Cold side: regen coolant on the outside of the block (same gap, larger radius so slower flow)
A_flow_out = 2 * np.pi * r_block_out * gap
v_out = mdot_fuel / (rho_cool * A_flow_out)
Re_out = rho_cool * v_out * D_h / mu_cool
h_c_out = 0.023 * Re_out**0.8 * Pr_cool**0.4 * k_cool / D_h

# Split the ring into N layers (finite volumes), per unit length of engine
N = 40
r_faces = np.linspace(R_throat, r_block_out, N + 1) # layer boundaries (m)
r_centers = 0.5 * (r_faces[:-1] + r_faces[1:]) # layer centers (m)
V_layer = np.pi * (r_faces[1:]**2 - r_faces[:-1]**2) # layer volume per meter of length (m^2)
G = k_copper * 2 * np.pi * r_faces[1:-1] / (r_centers[1:] - r_centers[:-1]) # conductance between neighbouring layers (W/(m*K))
dt = 0.25 * (r_faces[1] - r_faces[0])**2 / alpha_cu # time step small enough to stay stable (s)

T_cu = np.full(N, 298.0) # copper starts at room temperature (K)
t = 0.0
t_limit = None # time the hot surface reaches the copper limit
time_hist, T_surface_hist, T_outer_hist = [], [], []
for n in range(int(burn_time / dt)):
    heat = np.zeros(N) # net heat into each layer (W per meter)
    flow = G * (T_cu[:-1] - T_cu[1:]) # heat flowing outward between layers
    heat[:-1] -= flow
    heat[1:] += flow
    heat[0] += h_g_throat * (T_seen_throat - T_cu[0]) * 2 * np.pi * R_throat # hot gas into inner surface
    heat[-1] -= h_c_out * (T_cu[-1] - T_coolant) * 2 * np.pi * r_block_out # coolant out of outer surface
    T_cu += heat * dt / (rho_cu * cp_cu * V_layer)
    t += dt
    if t_limit is None and T_cu[0] >= T_copper_limit:
        t_limit = t
    if n % 200 == 0:
        time_hist.append(t); T_surface_hist.append(T_cu[0]); T_outer_hist.append(T_cu[-1])

print("\n-- Copper Heat-Sink Throat --")
print(f"Block:                       copper ring, r = {R_throat*1e3:.2f} to {r_block_out*1e3:.1f} mm")
print(f"Burn Simulated:              {burn_time:.1f} s at full thrust")
print(f"Hot Surface at End:          {T_cu[0]:.0f} K (limit {T_copper_limit} K)")
print(f"Outer Surface at End:        {T_cu[-1]:.0f} K")
print(f"Reaches Limit At:            {f'{t_limit:.1f} s' if t_limit else 'never (survives the burn)'}")

# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 6))
# plt.plot(time_hist, T_surface_hist, label='Hot (inner) surface')
# plt.plot(time_hist, T_outer_hist, label='Outer surface')
# plt.axhline(T_copper_limit, linestyle='--', color='red', label='Copper limit')
# plt.xlabel('Time (s)')
# plt.ylabel('Temperature (K)')
# plt.title('Copper throat heat soak')
# plt.grid(True)
# plt.legend()
# plt.show()
