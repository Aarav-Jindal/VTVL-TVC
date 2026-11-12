import numpy as np
from scipy.optimize import brentq
from CoolProp.CoolProp import PropsSI
import design as d

# Inputs (from design_values.json; geometry and gas properties come from cea_analysis.py)
psiToPa = 6894.76
Pc_Pa = d.get("Pc_psi") * psiToPa # chamber pressure (Pa)
mdot_fuel = d.get("mdot_fuel_kgps") # kg/s
cstar_real = d.get("cstar_real_mps") # m/s
mu_chosen = d.get("mu_gas_Pa_s") # Pa·s
cp_chosen = d.get("cp_gas_J_kgK") # J/(kg·K)
Pr_chosen = d.get("Pr_gas")
gamma_chosen = d.get("gamma")
T_0 = d.get("Tc_K") # chamber temperature at chosen mixture ratio (K)
R_chamber = d.get("R_chamber_mm") / 1e3 # m
R_throat = d.get("R_throat_mm") / 1e3 # m
R_exit = d.get("R_exit_mm") / 1e3 # m
L_cyl = d.get("L_cyl_mm") / 1e3 # m
L_converging = d.get("L_converging_mm") / 1e3 # m
L_nozzle = d.get("L_nozzle_mm") / 1e3 # m
L_total = d.get("L_total_mm") / 1e3 # m
convering_angle = d.get("converging_angle_deg") # degrees
chosen_angle = d.get("diverging_angle_deg") # degrees
Dt = 2 * R_throat # m

# Bartz at the throat
T_wg = d.get("T_wg_K") # hot-side wall temperature at the throat (K), placeholder
M = 1 # mach number (always 1 at the throat)
R_curv = 1.5 * R_throat # throat radius of curvature = R_up (m)
area_ratio = 1 # At/A (1 at the throat)
correctionFactor = 1/((0.5*(T_wg/T_0)*(1+(gamma_chosen-1)/2 * M**2)+1/2)**0.68*(1+(gamma_chosen-1)/2 * M**2)**0.12)
h_g = (0.026/Dt**0.2) * (mu_chosen**0.2 * cp_chosen / Pr_chosen**0.6) * (Pc_Pa/cstar_real)**0.8 * (Dt/R_curv)**0.1 * area_ratio**0.9 * correctionFactor # W/(m^2·K)
r = Pr_chosen**(1/3)
T_aw = T_0 * (1 + r * (gamma_chosen-1)/2 * M**2)/(1 + (gamma_chosen-1)/2 * M**2) # adiabatic wall temperature (K)
q = h_g * (T_aw - T_wg) # heat flux (W/m^2)

print("\n-- Heat Transfer --")
print(f"Adiabatic Wall Temperature:  {T_aw:.2f} K")
print(f"Heat Flux:                   {q/1e6:.2f} MW/m^2")

# Heat Flux along engine calculations

# Step 1: radius of the inner wall at any x (m), x measured from the injector face
# Same shape as the CAD: chamber, 45 deg cone, throat fillet R_up, short flat at Rt, fillet R_down, 15 deg cone.
# The two fillets need a flat of length L_throat between the two "virtual" sharp corners, which pushes the
# nozzle downstream by L_throat compared with a single sharp throat (so the total length grows by L_throat).
R_up = 1.5 * R_throat # upstream throat fillet radius (m)
R_down = 0.4 * R_throat # downstream throat fillet radius (m)
a_conv = np.radians(convering_angle)
a_div = np.radians(chosen_angle)
t_up = R_up * np.tan(a_conv / 2) # tangent length of the upstream fillet (m)
t_down = R_down * np.tan(a_div / 2) # tangent length of the downstream fillet (m)
L_throat = t_up + t_down + 0.01e-3 # flat between the virtual sharp corners, same equation as the CAD (m)
x_A = L_cyl + L_converging # first virtual sharp corner, end of the 45 deg cone (m)
x_B = x_A + L_throat # second virtual sharp corner, start of the 15 deg cone (m)
L_total = x_B + L_nozzle # injector face to nozzle exit, including the throat flat (m)
x_up_start = x_A - t_up * np.cos(a_conv) # upstream fillet starts on the cone (m)
x_up_end = x_A + t_up # upstream fillet meets the flat (m)
x_down_start = x_B - t_down # downstream fillet leaves the flat (m)
x_down_end = x_B + t_down * np.cos(a_div) # downstream fillet meets the 15 deg cone (m)

def radius(x):
    if x <= L_cyl:
        return R_chamber
    elif x <= x_up_start:
        return R_chamber - (x - L_cyl) * np.tan(a_conv)
    elif x <= x_up_end:
        return R_throat + R_up - np.sqrt(max(R_up**2 - (x - x_up_end)**2, 0))
    elif x <= x_down_start:
        return R_throat
    elif x <= x_down_end:
        return R_throat + R_down - np.sqrt(max(R_down**2 - (x - x_down_start)**2, 0))
    else:
        return R_throat + (x - x_B) * np.tan(a_div)

def slope(x):
    # dr/dx of the wall
    h = 1e-7
    return (radius(x + h) - radius(x - h)) / (2 * h)

# Wall angle at any x (degrees), needed to get the slice length along the wall
def wall_angle(x):
    return np.degrees(np.arctan(abs(slope(x))))

# Step 2: stations along the engine
dx = 5e-4 # 0.5 mm spacing (m)
x_stations = np.arange(0, L_total, dx)
x_throat = 0.5 * (x_up_end + x_down_start) # middle of the throat flat (m)

# Step 3: Mach number from area ratio A/At
def area_mach(M):
    return (1/M) * ((2/(gamma_chosen+1)) * (1 + (gamma_chosen-1)/2 * M**2))**((gamma_chosen+1)/(2*(gamma_chosen-1)))

def mach(areaRatio, supersonic):
    def f(M): return area_mach(M) - areaRatio
    if supersonic:
        return brentq(f, 1.0, 10.0)
    else:
        return brentq(f, 1e-6, 1.0)

# Step 4: Bartz heat flux at every station
r_recovery = Pr_chosen**(1/3) # recovery factor
r_list, M_list, q_list = [], [], []
h_g_list, T_aw_list = [], [] # stored for the film cooling model
Q_total = 0 # total heat into the wall (W)
for each_x in x_stations:
    r_wall = radius(each_x)
    A_over_At = (r_wall / R_throat)**2
    M_x = mach(A_over_At, each_x > x_throat)
    k = 1 + (gamma_chosen-1)/2 * M_x**2
    sigma_x = 1/((0.5*(T_wg/T_0)*k + 0.5)**0.68 * k**0.12)
    h_g_x = (0.026/Dt**0.2) * (mu_chosen**0.2 * cp_chosen / Pr_chosen**0.6) * (Pc_Pa/cstar_real)**0.8 * (Dt/R_curv)**0.1 * (1/A_over_At)**0.9 * sigma_x
    T_aw_x = T_0 * (1 + r_recovery*(gamma_chosen-1)/2 * M_x**2) / k
    q_x = h_g_x * (T_aw_x - T_wg) # W/m^2

    # Step 5: add this slice's heat to the total
    ds = dx / np.cos(np.radians(wall_angle(each_x))) # length along the wall (m)
    Q_total += q_x * 2 * np.pi * r_wall * ds

    r_list.append(r_wall); M_list.append(M_x); q_list.append(q_x)
    h_g_list.append(h_g_x); T_aw_list.append(T_aw_x)

q_list = np.array(q_list)
i_max = np.argmax(q_list)

# Step 6: can the ethanol absorb it?
cp_ethanol = 2500 # J/(kg·K), placeholder
dT_coolant = 150 - 25 # allowed ethanol temperature rise (K), placeholder
coolantCapacity = mdot_fuel * cp_ethanol * dT_coolant # W

print("\n-- Regen Heat Load --")
print(f"Peak heat flux:              {q_list[i_max]/1e6:.2f} MW/m^2 at x = {x_stations[i_max]*1e3:.1f} mm (throat at {x_throat*1e3:.1f} mm)")
print(f"Total heat into wall:        {Q_total/1e3:.2f} kW")
print(f"Ethanol can absorb:          {coolantCapacity/1e3:.2f} kW")
print(f"Enough cooling?              {'YES' if Q_total < coolantCapacity else 'NO'}")

# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 6))
# plt.plot(x_stations*1e3, q_list/1e6)
# plt.axvline(x_throat*1e3, linestyle='--', color='gray', label='Throat')
# plt.xlabel('Distance from injector face (mm)')
# plt.ylabel('Heat flux (MW/m^2)')
# plt.title('Bartz heat flux along the engine')
# plt.grid(True)
# plt.legend()
# plt.show()

# Film cooling Calculation
# Ethanol boiling properties at chamber pressure (Q = vapor quality: 0 = all liquid, 1 = all vapor)
T_boil = PropsSI('T', 'P', Pc_Pa, 'Q', 0, 'Ethanol') # boiling temperature (K)
h_liq = PropsSI('H', 'P', Pc_Pa, 'Q', 0, 'Ethanol') # enthalpy of liquid just at boiling (J/kg)
h_vap = PropsSI('H', 'P', Pc_Pa, 'Q', 1, 'Ethanol') # enthalpy of vapor just after boiling (J/kg)
h_fg = h_vap - h_liq # heat of vaporization (J/kg)

T_film_in = 423.15 # film ethanol enters after the cooling jacket, at 150 C (K)
h_in = PropsSI('H', 'T', T_film_in, 'P', Pc_Pa, 'Ethanol') # enthalpy of film entering (J/kg)
E = h_vap - h_in # heat 1 kg of film absorbs before it is fully boiled (J/kg)

print("\n-- Film Cooling (Ethanol at Chamber Pressure) --")
print(f"Boiling Temperature:         {T_boil:.1f} K")
print(f"Heat of Vaporization:        {h_fg/1e3:.1f} kJ/kg")
print(f"Film Heat Capacity (E):      {E/1e3:.1f} kJ/kg")

# Step 2: effective film flow (only part of the film stays on the wall; the rest is torn off into the core)
filmFraction = d.get("film_fraction") # fraction of the ethanol sent to the film, f
eta_film = d.get("film_eta") # fraction of the film that actually stays on the wall (placeholder, typical value)
mdot_film = filmFraction * mdot_fuel # film flow (kg/s)
mdot_film_eff = eta_film * mdot_film # film flow that actually protects the wall (kg/s)
filmHeatBudget = mdot_film_eff * E # heat the film absorbs per second before it is gone (W)

print(f"Film Fraction (f):           {filmFraction*100:.0f} % of fuel")
print(f"Film Flow:                   {mdot_film*1e3:.2f} g/s")
print(f"Effective Film Flow:         {mdot_film_eff*1e3:.2f} g/s")
print(f"Film Heat Budget:            {filmHeatBudget/1e3:.2f} kW")

# Step 3: how far down the engine the liquid film survives
# heat into each wall slice = q * (2*pi*r) * ds, using the stations from the regen loop
ds_list = np.array([dx / np.cos(np.radians(wall_angle(x))) for x in x_stations]) # slice length along the wall (m)
sliceHeat = q_list * 2 * np.pi * np.array(r_list) * ds_list # heat into each slice (W)
heatSoFar = np.cumsum(sliceHeat) # running total of heat from the injector face (W)

if heatSoFar[-1] <= filmHeatBudget: # film never runs out
    i_film_end = len(x_stations) - 1
else:
    i_film_end = np.argmax(heatSoFar >= filmHeatBudget) # first station where the film is used up
x_film_end = x_stations[i_film_end]

print(f"Liquid Film Ends At:         x = {x_film_end*1e3:.1f} mm (chamber ends at {L_cyl*1e3:.1f} mm, throat at {x_throat*1e3:.1f} mm)")
print(f"Heat Covered by Film:        {heatSoFar[i_film_end]/1e3:.2f} kW of {Q_total/1e3:.2f} kW")

# Vapor film: after the liquid boils off, the ethanol vapor stays near the wall as a cooler layer.
# The wall then "sees" the film temperature instead of the hot core gas (T_aw).
# The core gas heats the film as it travels: m*cp*dT_f = h_g*(T_aw - T_f)*dA, so the film warms toward T_aw
# and its protection fades with distance.
# First-pass model: no extra mixing beyond eta_film, film heating counted from the core side only (conservative).
cp_film_vapor = d.get("cp_film_vapor") # ethanol vapor cp (J/(kg*K)), placeholder average for ~500-1000 K
T_film = T_boil # film temperature where the liquid ends (K)
q_wall_film = np.zeros(len(x_stations)) # wall heat flux with film cooling (W/m^2)
for j in range(len(x_stations)):
    if j <= i_film_end:
        continue # wall under liquid film: its heat goes into boiling the film, not into the regen
    T_seen = min(T_film, T_aw_list[j]) # gas temperature the wall sees
    q_wall_film[j] = max(h_g_list[j] * (T_seen - T_wg), 0) # no heat flow into the wall while film is cooler than it
    dA = 2 * np.pi * r_list[j] * ds_list[j] # wall area of this slice (m^2)
    T_film += h_g_list[j] * (T_aw_list[j] - T_film) * dA / (mdot_film_eff * cp_film_vapor) # film heats up

Q_wall_film = np.sum(q_wall_film * 2 * np.pi * np.array(r_list) * ds_list) # total wall heat with film (W)
i_throat_station = np.argmin(np.abs(x_stations - x_throat))

print(f"Wall Heat Flux at Throat:    {q_wall_film[i_throat_station]/1e6:.2f} MW/m^2 (was {q_list[i_throat_station]/1e6:.2f} without film)")
print(f"Total Wall Heat with Film:   {Q_wall_film/1e3:.2f} kW (was {Q_total/1e3:.2f} kW without film)")

# Regen channel design
# N printed channels of constant width around a copper liner, running the full length of the engine.
# Ethanol enters at the nozzle exit and flows toward the injector (counterflow).
# Where straight channels can't keep the wall cool enough, the channels twist into a helix
# (helix angle beta, measured from the engine axis):
#   - neighbouring channels sit closer together (rib across the flow = 2*pi*r*cos(beta)/N - w),
#     so more of the hot wall touches coolant
#   - each channel gets longer (ds / cos(beta)), so the pressure drop goes up
# Step 1: at each station, the smallest twist that can keep the wall at T_wg_max with a printable channel.
# Step 2: smooth the twist so it changes gradually along the wall (no kinks in the channels).
# Step 3: with that twist, the deepest channel (lowest pressure drop) that still keeps the wall at T_wg_max.
k_copper = d.get("k_liner") # liner thermal conductivity (W/(m*K)), CuCrZr
P_tank_psi = d.get("P_tank_psi") # tank pressure (psia)
P_jacket = (P_tank_psi - d.get("feed_losses_psi")) * psiToPa # ethanol pressure entering the jacket (Pa): tank - lines/valves
N_ch = int(d.get("n_channels")) # number of channels
w_ch = d.get("channel_width_m") # channel width (m), constant along the engine
t_hot = d.get("hot_wall_m") # hot wall thickness between gas and channel bottom (m), placeholder
T_wg_max = T_wg # design hot-side wall temperature (K), same as Bartz assumed
T_wall_max = d.get("T_wall_max_K") # most the wall may reach, even with print errors (K)
roughness = d.get("channel_roughness_m") # as-printed channel surface roughness (m), placeholder
depth_min = d.get("channel_depth_min_m") # shallowest printable channel (m)
rib_min = d.get("rib_min_m") # thinnest printable rib (m)
beta_max = d.get("helix_angle_max_deg") # steepest twist allowed (deg)
beta_ramp = d.get("helix_ramp_deg_per_mm") # fastest the twist may change along the wall (deg per mm)
overhang_depth = d.get("helix_overhang_depth_max_m") # deepest channel allowed above 45 deg of twist (m)
print_tol = d.get("print_tolerance_m") # print size error for the wall temperature check (m)
depths = np.arange(depth_min, 3.0e-3, 0.01e-3) # candidate channel depths (m)
betas = np.arange(0, beta_max + 0.5, 1.0) # candidate helix angles (deg)
T_cool_in = d.get("T_coolant_in_K") # ethanol temperature entering the jacket (K)
P_crit_eth = PropsSI('pcrit', 'Ethanol') # ethanol critical pressure (Pa)
if P_jacket < P_crit_eth:
    T_sat_jacket = PropsSI('T', 'P', P_jacket, 'Q', 0, 'Ethanol') # ethanol boiling temperature in the jacket (K)
else:
    T_sat_jacket = None # above the critical pressure ethanol cannot boil

def channel_h_eff(d, w, land, props):
    # Effective heat transfer coefficient at the channel's hot side, including the ribs acting as fins
    rho, mu, k_f, Pr_f = props
    A = w * d # channel flow area (m^2)
    Dh = 4 * A / (2 * (w + d)) # hydraulic diameter (m)
    v = (mdot_fuel / N_ch) / (rho * A) # coolant velocity (m/s)
    Re = rho * v * Dh / mu
    f = (-1.8 * np.log10((roughness / Dh / 3.7)**1.11 + 6.9 / Re))**-2 # Darcy friction factor (Haaland)
    if Re > 2300:
        Nu = (f/8) * (Re - 1000) * Pr_f / (1 + 12.7 * np.sqrt(f/8) * (Pr_f**(2/3) - 1)) # Gnielinski
        Nu = max(Nu, 4.36)
    else:
        Nu = 4.36 # laminar
        f = 64 / Re
    h_c = Nu * k_f / Dh
    m_fin = np.sqrt(2 * h_c / (k_copper * land))
    eta_fin = np.tanh(m_fin * d) / (m_fin * d) # rib fin efficiency
    h_eff = h_c * (w + 2 * eta_fin * d) / (w + land) # per unit area of channel pitch
    return h_eff, v, Dh, Re, f

def rib_width(r_cold, beta):
    # rib between neighbouring channels, measured across the flow (m)
    return 2 * np.pi * r_cold * np.cos(np.radians(beta)) / N_ch - w_ch

def printable_depths(beta):
    # above 45 deg of twist the channel roof faces down while printing: keep it shallow enough to bridge
    return depths if beta <= 45 else depths[depths <= overhang_depth]

n_st = len(x_stations)
s_mm = (np.cumsum(ds_list) - ds_list / 2) * 1e3 # distance along the wall from the injector face (mm)

# Coolant temperature and properties along the engine. They only depend on the heat picked up,
# not on the channel shape, so they're found first.
props_list = [None] * n_st; cp_list = np.zeros(n_st)
T_cool_in_slice = np.zeros(n_st); T_cool_ch = np.zeros(n_st)
T_cool_now = T_cool_in
for j in reversed(range(n_st)): # march from nozzle exit to injector
    props_list[j] = (PropsSI('D', 'T', T_cool_now, 'P', P_jacket, 'Ethanol'),
                     PropsSI('V', 'T', T_cool_now, 'P', P_jacket, 'Ethanol'),
                     PropsSI('L', 'T', T_cool_now, 'P', P_jacket, 'Ethanol'),
                     PropsSI('PRANDTL', 'T', T_cool_now, 'P', P_jacket, 'Ethanol'))
    cp_list[j] = PropsSI('C', 'T', T_cool_now, 'P', P_jacket, 'Ethanol')
    T_cool_in_slice[j] = T_cool_now
    T_cool_now += q_wall_film[j] * 2 * np.pi * r_list[j] * ds_list[j] / (mdot_fuel * cp_list[j]) # coolant heats up
    T_cool_ch[j] = T_cool_now

def h_needed(j):
    # h_eff that keeps the hot wall at T_wg_max (heat per length is the same through every layer)
    q_j = q_wall_film[j]; r_hot = r_list[j]
    T_wc_allowed = T_wg_max - q_j * t_hot / k_copper # coolant-side wall temperature allowed (K)
    if T_wc_allowed <= T_cool_in_slice[j]:
        return np.inf # impossible even with perfect cooling
    return q_j * (r_hot / (r_hot + t_hot)) / (T_wc_allowed - T_cool_in_slice[j])

# Step 1: smallest twist that works at each station
beta_need = np.zeros(n_st); beta_cap = np.zeros(n_st)
for j in range(n_st):
    r_cold = r_list[j] + t_hot
    cos_cap = min(N_ch * (w_ch + rib_min) / (2 * np.pi * r_cold), 1.0) # twist where the rib gets too thin
    beta_cap[j] = min(beta_max, np.degrees(np.arccos(cos_cap)))
    if q_wall_film[j] <= 0:
        continue # no heat (liquid film region): straight
    h_need = h_needed(j)
    beta_need[j] = beta_cap[j] # if nothing works, twist as far as the ribs allow
    for beta in betas[betas <= beta_cap[j]]:
        land = rib_width(r_cold, beta)
        if any(channel_h_eff(dd, w_ch, land, props_list[j])[0] >= h_need for dd in printable_depths(beta)):
            beta_need[j] = beta
            break

# Step 2: smooth it. Each station takes the most twist any neighbour needs, less beta_ramp per mm of distance.
beta_ch = np.array([np.max(beta_need - beta_ramp * np.abs(s_mm - s_mm[j])) for j in range(n_st)])
beta_ch = np.minimum(beta_ch, beta_cap)

# Step 3: deepest printable channel that keeps the wall at T_wg_max with that twist
depth_ch = np.zeros(n_st); land_ch = np.zeros(n_st); vel_ch = np.zeros(n_st)
T_wall_ch = np.zeros(n_st); T_wall_tol = np.zeros(n_st)
dP_total = 0.0
infeasible = []
for j in range(n_st):
    r_hot = r_list[j]
    r_cold = r_hot + t_hot # radius at channel bottom (m)
    land_j = rib_width(r_cold, beta_ch[j])
    props = props_list[j]
    q_j = q_wall_film[j]
    allowed = printable_depths(beta_ch[j])
    if q_j <= 0:
        d_ch = allowed[-1] # no heat (liquid film region): deepest channel
    else:
        h_need = h_needed(j)
        ok = [dd for dd in allowed if channel_h_eff(dd, w_ch, land_j, props)[0] >= h_need]
        if ok:
            d_ch = max(ok) # deepest channel that is still good enough
        else:
            d_ch = allowed[0]; infeasible.append(j) # even the shallowest channel isn't enough

    h_eff, v, Dh, Re, f = channel_h_eff(d_ch, w_ch, land_j, props)
    T_wall_ch[j] = T_cool_in_slice[j] + q_j * (r_hot / r_cold) / h_eff + q_j * t_hot / k_copper # actual hot-side wall temperature (K)
    # Same channel printed print_tol too deep and too wide (slower coolant, hotter wall)
    h_tol = channel_h_eff(d_ch + print_tol, w_ch + print_tol, land_j - print_tol, props)[0]
    T_wall_tol[j] = T_cool_in_slice[j] + q_j * (r_hot / r_cold) / h_tol + q_j * t_hot / k_copper
    dP_total += f * (ds_list[j] / np.cos(np.radians(beta_ch[j])) / Dh) * props[0] * v**2 / 2 # pressure drop in this slice (Pa)
    depth_ch[j] = d_ch; land_ch[j] = land_j; vel_ch[j] = v
T_cool_now = T_cool_ch[0] # coolant leaving the jacket at the injector end (K)

# Where the channel sits around the engine: it turns by tan(beta) * ds / r at its mid-depth radius
r_mid = np.array(r_list) + t_hot + depth_ch / 2
theta_ch = np.cumsum(np.tan(np.radians(beta_ch)) * ds_list / r_mid) # radians, 0 at the injector face
L_channel = np.sum(ds_list / np.cos(np.radians(beta_ch))) # length of one channel (m)
i_helix_start = int(np.argmax(beta_ch > 0)) if beta_ch.max() > 0 else n_st - 1

P_inj_in_eth = P_jacket/psiToPa - dP_total/psiToPa # ethanol pressure at injector after the jacket (psia)
dP_inj_eth = P_inj_in_eth - Pc_Pa/psiToPa # what's left for the ethanol injector (psi)

print("\n-- Regen Channel Design --")
print(f"Channels:                    {N_ch}, width {w_ch*1e3:.2f} mm, hot wall {t_hot*1e3:.2f} mm")
print(f"Rib at Throat / Chamber:     {land_ch[i_throat_station]*1e3:.2f} / {land_ch[0]*1e3:.2f} mm")
print(f"Depth at Throat (min/max):   {depth_ch[i_throat_station]*1e3:.2f} mm ({depth_ch.min()*1e3:.2f} / {depth_ch.max()*1e3:.2f} mm)")
print(f"Velocity at Throat (max):    {vel_ch[i_throat_station]:.1f} m/s ({vel_ch.max():.1f} m/s)")
print(f"Helix Angle:                 straight up to x = {x_stations[i_helix_start]*1e3:.1f} mm, {beta_ch[i_throat_station]:.0f} deg at throat, {beta_ch.max():.0f} deg max")
print(f"Channel Twist / Length:      {np.degrees(theta_ch[-1]):.0f} deg around the engine, {L_channel*1e3:.0f} mm long")
print(f"Max Hot Wall Temperature:    {T_wall_ch.max():.0f} K (design {T_wg_max} K)")
print(f"Wall If Printed Oversize:    {T_wall_tol.max():.0f} K with channels {print_tol*1e3:.2f} mm deeper and wider (max allowed {T_wall_max} K)")
print(f"Coolant Out (to injector):   {T_cool_now:.0f} K " + (f"(boils at {T_sat_jacket:.0f} K in jacket)" if T_sat_jacket else f"(jacket above critical pressure {P_crit_eth/psiToPa:.0f} psi: no boiling)"))
print(f"Jacket Pressure Drop:        {dP_total/psiToPa:.1f} psi")
print(f"Ethanol Injector Drop Left:  {dP_inj_eth:.1f} psi = {dP_inj_eth/(Pc_Pa/psiToPa)*100:.0f} % of Pc (need >= ~20 %)")
print(f"Infeasible Stations:         {len(infeasible)}" + (f" (x = {x_stations[min(infeasible)]*1e3:.1f} to {x_stations[max(infeasible)]*1e3:.1f} mm)" if infeasible else ""))

# Save results for the other scripts
d.save_many({
    "q_throat_nofilm_W_m2": (q_list[i_throat_station], "W/m^2"),
    "q_throat_film_W_m2": (q_wall_film[i_throat_station], "W/m^2"),
    "h_g_throat_W_m2K": (h_g_list[i_throat_station], "W/(m^2*K)"),
    "Q_wall_nofilm_W": (Q_total, "W"),
    "Q_wall_film_W": (Q_wall_film, "W"),
    "mdot_film_kgps": (mdot_film, "kg/s"),
    "x_liquid_film_end_mm": (x_film_end*1e3, "mm"),
    "T_coolant_out_K": (T_cool_now, "K"),
    "P_jacket_in_psi": (P_jacket/psiToPa, "psia"),
    "jacket_dP_psi": (dP_total/psiToPa, "psi"),
    "P_ethanol_injector_in_psi": (P_inj_in_eth, "psia"),
    "dP_ethanol_injector_psi": (dP_inj_eth, "psi"),
    "max_wall_T_K": (T_wall_ch.max(), "K"),
    "channel_depth_throat_mm": (depth_ch[i_throat_station]*1e3, "mm"),
    "channel_depth_min_mm": (depth_ch.min()*1e3, "mm"),
    "channel_depth_max_mm": (depth_ch.max()*1e3, "mm"),
    "regen_infeasible_stations": (len(infeasible), "-"),
    "max_wall_T_print_error_K": (T_wall_tol.max(), "K"),
    "helix_angle_throat_deg": (beta_ch[i_throat_station], "deg"),
    "helix_angle_peak_deg": (beta_ch.max(), "deg"),
    "x_helix_start_mm": (x_stations[i_helix_start]*1e3, "mm"),
    "channel_twist_deg": (np.degrees(theta_ch[-1]), "deg"),
    "channel_length_mm": (L_channel*1e3, "mm"),
    "L_gas_total_mm": (L_total*1e3, "mm"),
    "L_throat_flat_mm": (L_throat*1e3, "mm"),
    "x_throat_mm": (x_throat*1e3, "mm"),
    "R_up_mm": (R_up*1e3, "mm"),
    "R_down_mm": (R_down*1e3, "mm"),
}, "heat_transfer.py")

# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 6))
# plt.plot(x_stations*1e3, depth_ch*1e3, label='Channel depth (mm)')
# plt.plot(x_stations*1e3, land_ch*1e3, label='Rib width (mm)')
# plt.axvline(x_throat*1e3, linestyle='--', color='gray', label='Throat')
# plt.xlabel('Distance from injector face (mm)')
# plt.ylabel('mm')
# plt.title('Regen channel geometry')
# plt.grid(True)
# plt.legend()
# plt.show()

# Export curves for SolidWorks (Insert > Curve > Curve Through XYZ Points)
# Each file is "x y z" in mm: x along the engine axis from the injector face, y = radius, z = 0.
# Offsets are taken perpendicular to the wall (not straight up), so the 0.8 mm hot wall is 0.8 mm
# thick on the cones too. Throat fillets are not in these profiles: for the channel floor, use
# Offset Entities (hot wall thickness) on your filleted gas-path sketch instead.
import os

def wall_normal(x):
    # outward unit normal of the gas-path wall at x (axial, radial), from the wall slope
    s_ = slope(x)
    n = np.hypot(1.0, s_)
    return -s_ / n, 1.0 / n

# Jacket wall: the metal outside the channels holds the coolant pressure (hoop stress, like the tanks)
# t = P * SF * r / sigma, worst case at the biggest radius (the chamber), rounded up to 0.1 mm
sigma_jacket = d.get("jacket_yield_Pa") # Pa
SF_jacket = d.get("jacket_safety_factor")
r_jacket_inner = np.array([r_list[j] + t_hot + depth_ch[j] for j in range(len(x_stations))]) # rib tops = jacket inside (m)
t_jacket_required = P_jacket * SF_jacket * r_jacket_inner / sigma_jacket # required wall at each station (m)
t_jacket = max(np.ceil(t_jacket_required.max() * 1e4) / 1e4, d.get("jacket_min_m")) # one constant wall (m)
stock = d.get("machining_stock_m") # extra metal on the gas side, machined off after printing (m)

print("\n-- Jacket Wall --")
print(f"Coolant Pressure:            {P_jacket/psiToPa:.0f} psi")
print(f"Required (worst, chamber):   {t_jacket_required.max()*1e3:.2f} mm at r = {r_jacket_inner.max()*1e3:.2f} mm")
print(f"Chosen Jacket Wall:          {t_jacket*1e3:.2f} mm (min {d.get('jacket_min_m')*1e3:.1f} mm)")
print(f"Machining Stock (gas side):  {stock*1e3:.2f} mm")

d.save_many({
    "t_jacket_mm": (t_jacket*1e3, "mm"),
    "jacket_outer_r_chamber_mm": ((r_jacket_inner[0] + t_jacket)*1e3, "mm"),
}, "heat_transfer.py")

curve_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CAD", "curves")
os.makedirs(curve_folder, exist_ok=True)
curves = {"gas_path.txt": [], "gas_path_printed.txt": [], "channel_floor.txt": [], "liner_outer.txt": [], "jacket_outer.txt": []}
for j, x in enumerate(x_stations):
    nx, ny = wall_normal(x)
    r = r_list[j]
    curves["gas_path.txt"].append((x, r))
    curves["channel_floor.txt"].append((x + nx * t_hot, r + ny * t_hot))
    curves["liner_outer.txt"].append((x + nx * (t_hot + depth_ch[j]), r + ny * (t_hot + depth_ch[j])))
    curves["gas_path_printed.txt"].append((x - nx * stock, r - ny * stock)) # gas side with machining stock (inward)
    curves["jacket_outer.txt"].append((x + nx * (t_hot + depth_ch[j] + t_jacket), r + ny * (t_hot + depth_ch[j] + t_jacket)))
for name, pts in curves.items():
    with open(os.path.join(curve_folder, name), "w") as f:
        for x, y in pts:
            f.write(f"{x*1e3:.4f}\t{y*1e3:.4f}\t0\n")

# 3D helix of one channel: its floor centre line and rib-top centre line, turned to theta around the axis.
# Clock 0 = +y, turning toward +z. Sweep the channel along the floor line, guided by the top line.
helix = {"channel_floor_helix.txt": t_hot, "channel_top_helix.txt": None}
for name, off in helix.items():
    with open(os.path.join(curve_folder, name), "w") as f:
        for j, x in enumerate(x_stations):
            nx, ny = wall_normal(x)
            o = t_hot if off is not None else t_hot + depth_ch[j]
            xx, rr = x + nx * o, r_list[j] + ny * o
            f.write(f"{xx*1e3:.4f}\t{rr*np.cos(theta_ch[j])*1e3:.4f}\t{rr*np.sin(theta_ch[j])*1e3:.4f}\n")

print("\n-- CAD Curves Exported --")
print(f"Folder:                      {curve_folder}")
print(f"Files:                       {', '.join(list(curves) + list(helix))} ({len(x_stations)} points each, mm)")
print(f"Liner material:              {d.get('liner_material')}")
