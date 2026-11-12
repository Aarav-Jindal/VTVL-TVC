# Printed chamber hardware: coolant manifolds, flange, seals and bolts.
# Checks the sizes in design_values.json and prints every dimension needed for the CAD.
# x is measured from the injector face (flange face = x 0), r from the engine axis.
# Run heat_transfer.py first: this reads its exported channel curves and pressures.
import os
import numpy as np
from CoolProp.CoolProp import PropsSI
import design as d

psiToPa = 6894.76
mdot_fuel = d.get("mdot_fuel_kgps") # all the ethanol goes through the jacket (kg/s)
Pc_Pa = d.get("Pc_psi") * psiToPa
P_jacket_in = d.get("P_jacket_in_psi") * psiToPa # coolant pressure at the inlet manifold (Pa)
P_jacket_out = P_jacket_in - d.get("jacket_dP_psi") * psiToPa # coolant pressure at the outlet manifold (Pa)
T_in = d.get("T_coolant_in_K")
T_out = d.get("T_coolant_out_K")
rho_in = PropsSI('D', 'T', T_in, 'P', P_jacket_in, 'Ethanol') # cold ethanol (kg/m^3)
rho_out = PropsSI('D', 'T', T_out, 'P', P_jacket_out, 'Ethanol') # hot ethanol (kg/m^3)
n_ch = int(d.get("n_channels"))
w_ch = d.get("channel_width_m")
L_total = d.get("L_gas_total_mm") / 1e3 # injector face to nozzle exit, including the throat flat (m)
R_chamber = d.get("R_chamber_mm") / 1e3

# Channel floor and rib-top radius along the engine, from heat_transfer.py's exported curves (mm -> m)
curve_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CAD", "curves")
floor = np.loadtxt(os.path.join(curve_folder, "channel_floor.txt"))[:, :2] / 1e3
ribtop = np.loadtxt(os.path.join(curve_folder, "liner_outer.txt"))[:, :2] / 1e3
jacket = np.loadtxt(os.path.join(curve_folder, "jacket_outer.txt"))[:, :2] / 1e3
def floor_r(x): return np.interp(x, floor[:, 0], floor[:, 1])
def ribtop_r(x): return np.interp(x, ribtop[:, 0], ribtop[:, 1])
def jacket_r(x): return np.interp(x, jacket[:, 0], jacket[:, 1])

# Rule for even flow: velocity in a manifold should be well below the velocity in the channels it feeds,
# so every channel sees nearly the same pressure.
# Both manifolds are rings sitting on top of the ribs: the channels run underneath them and open into them
# through their tops, so the hot wall is channel-cooled everywhere (slow manifold flow would barely cool it).
# The chamber prints flange-down, so each manifold roof slopes at least roof_angle from flat (no supports inside).
roof_tan = np.tan(np.radians(d.get("manifold_roof_angle_deg")))
rib_cut = d.get("manifold_rib_cut_m") # manifold inner surface sits this far into the rib tops (m)
def polygon_area(pts):
    x, y = np.array(pts).T
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

# Inlet manifold (nozzle end): cold ethanol comes in through one port, splits both ways around the ring,
# then drops into the channels. Its top is where the channels end.
x_in_top = L_total - d.get("inlet_manifold_end_gap_m") # top of the manifold = end of the channels (m)
x_in_floor = x_in_top - d.get("inlet_manifold_len_m") # floor: the inlet passage comes up through it (m)
r_in_outer = d.get("inlet_manifold_outer_r_m")
r_in_rib_floor = ribtop_r(x_in_floor) - rib_cut # inner edge (on the rib tops) at the floor and at the top (m)
r_in_rib_top = ribtop_r(x_in_top) - rib_cut
x_in_roof_out = x_in_top - (r_in_outer - r_in_rib_top) * roof_tan # roof drops outward from the top to the outer wall (m)
if x_in_roof_out > x_in_floor:
    in_poly = [(x_in_floor, r_in_rib_floor), (x_in_top, r_in_rib_top), (x_in_roof_out, r_in_outer), (x_in_floor, r_in_outer)]
else: # roof reaches the floor before the outer wall
    in_poly = [(x_in_floor, r_in_rib_floor), (x_in_top, r_in_rib_top), (x_in_floor, r_in_rib_top + (x_in_top - x_in_floor) / roof_tan)]
r_in_floor_outer = in_poly[-1][1] # outer edge of the floor (m)
A_in_manifold = polygon_area(in_poly) # cross-section (m^2)
v_in_manifold = (mdot_fuel / 2) / (rho_in * A_in_manifold) # half the flow goes each way round (m/s)
A_ch_in = n_ch * w_ch * (ribtop_r(x_in_floor) - floor_r(x_in_floor)) # all channel area just below the manifold (m^2)
v_ch_in = mdot_fuel / (rho_in * A_ch_in)
A_port = np.pi / 4 * d.get("inlet_port_drill_d_m")**2
v_port = mdot_fuel / (rho_in * A_port)

# Outlet manifold (injector end): hot ethanol rises out of the channels, goes both ways round the ring to one outlet port,
# then leaves through an external line to the fuel throttle valve on the injector (the valve sits after the jacket,
# so the jacket stays near tank pressure at every thrust). Flat floor (prints facing up), roof rising inward at
# roof_angle from the outer wall to the rib tops.
x_out_floor = d.get("outlet_manifold_x_m")
r_out_outer = d.get("outlet_manifold_outer_r_m")
r_out_inner = ribtop_r(x_out_floor) - rib_cut # on the rib tops (m)
x_out_top = x_out_floor + (r_out_outer - r_out_inner) * roof_tan # where the roof meets the rib tops (m)
out_poly = [(x_out_floor, r_out_inner), (x_out_floor, r_out_outer), (x_out_top, r_out_inner)]
A_out_manifold = polygon_area(out_poly)
v_out_manifold = (mdot_fuel / 2) / (rho_out * A_out_manifold) # half the flow comes each way round to the port (m/s)
A_ch_out = n_ch * w_ch * (ribtop_r(x_out_top) - floor_r(x_out_top))
v_ch_out = mdot_fuel / (rho_out * A_ch_out)
# Even flow, the real test: the pressure change round a manifold (about one dynamic head) next to the jacket drop
dP_jacket = d.get("jacket_dP_psi") * psiToPa
head_in = rho_in * v_in_manifold**2 / 2 / dP_jacket
head_out = rho_out * v_out_manifold**2 / 2 / dP_jacket

# Outlet port: boss on the jacket just behind the flange, clock 0 (same side as the injector's ethanol inlet).
# Radial passage from the port down to r_pass, then axial along the engine into the outlet manifold through its roof.
x_op = d.get("outlet_port_x_m")
op_r = d.get("outlet_port_boss_d_m") / 2
r_op_top = d.get("outlet_port_boss_top_r_m")
r_op_bottom = jacket_r(x_op) - 1e-3 # boss starts 1 mm inside the jacket wall (m)
r_op_pass = d.get("outlet_passage_r_m"); d_op_pass = d.get("outlet_passage_d_m")
v_op_pass = mdot_fuel / (rho_out * np.pi / 4 * d_op_pass**2)
x_ch_start = x_out_floor + 0.5e-3 # channels start just above the outlet manifold floor (m)

# Seals: two face-seal O-ring grooves on the flange face
gw = d.get("oring_groove_width_m")
r_inner_seal = d.get("inner_groove_id_m") / 2 + gw # pressure pushes the ring outward: seals at the groove OD (m)
r_outer_seal = d.get("outer_groove_id_m") / 2 + gw

# Bolts: pressure tries to push the injector off the chamber
# Chamber gas acts inside the inner seal; if that seal leaks, gas reaches the outer (backup) seal.
F_sep = Pc_Pa * np.pi * r_outer_seal**2 # separation load, inner seal leaking (N)
n_bolts = int(d.get("n_bolts"))
proof_load = d.get("bolt_proof_Pa") * d.get("bolt_stress_area_m2") # N
preload = 0.75 * proof_load # usual preload target (N)
joint_SF = n_bolts * preload / F_sep

# Flange bending: bolt load at the bolt circle, seal reaction at the outer seal, ring bending
r_bc = d.get("bolt_circle_d_m") / 2
t_flange = d.get("flange_thickness_m")
arm = r_bc - r_outer_seal # lever arm (m)
m_ring = n_bolts * preload * arm / (2 * np.pi * r_bc) # bending moment per metre of circumference (N*m/m)
sigma_flange = 6 * m_ring / t_flange**2 # Pa
sigma_allow = d.get("jacket_yield_Pa") / d.get("jacket_safety_factor")
t_flange_needed = np.sqrt(6 * m_ring / sigma_allow)

# Geometry checks
metal_groove_to_manifold = x_out_floor - d.get("oring_groove_depth_m") # between groove bottom and outlet manifold (m)
# Seal lands on the flange face: solid metal between the bore, the two grooves and the bolt holes (m)
land_bore_inner = d.get("inner_groove_id_m") / 2 - R_chamber
land_inner_outer = d.get("outer_groove_id_m") / 2 - r_inner_seal
land_outer_bolts = r_bc - d.get("bolt_clearance_d_m") / 2 - r_outer_seal
# Outlet roof to the outside: closest point is the corner where the flange top meets the jacket
corner = np.array([t_flange, jacket_r(t_flange)])
roof_dir = np.array([x_out_top - x_out_floor, r_out_inner - r_out_outer]); roof_dir /= np.linalg.norm(roof_dir)
rel = corner - np.array([x_out_floor, r_out_outer])
metal_outlet_roof = abs(rel[0] * roof_dir[1] - rel[1] * roof_dir[0]) # distance from the corner to the roof line (m)
metal_inlet_top = L_total - x_in_top # inlet manifold top to the nozzle exit face (m)
metal_inlet_outer = d.get("inlet_collar_outer_r_m") - r_in_outer # inlet manifold to the collar outside (m)
# Every channel must open into both manifolds: the straight inner edge has to stay under the curved rib tops
xs_in = np.linspace(x_in_floor, x_in_top, 200)
edge_in = np.interp(xs_in, [x_in_floor, x_in_top], [r_in_rib_floor, r_in_rib_top])
inlet_opens_channels = np.all(edge_in < ribtop_r(xs_in))
r_boss_top = d.get("inlet_port_boss_top_r_m")
x_port = L_total - d.get("inlet_port_x_from_exit_m") # port boss centre (m)
x_collar_start = L_total - d.get("inlet_collar_len_m")
boss_r = d.get("inlet_port_boss_d_m") / 2
boss_fits = (x_port - boss_r >= x_collar_start) and (x_port + boss_r <= L_total)
r_pass = d.get("inlet_passage_r_m"); d_pass = d.get("inlet_passage_d_m")
pass_inner = r_pass - d_pass / 2 # inner edge of the internal passage (m)
xs_pass = np.linspace(x_port, x_in_floor, 50) # the axial passage runs from the port to the manifold floor
metal_under_passage = pass_inner - np.max(jacket_r(xs_pass)) # metal between jacket outer surface and passage, worst point along the run (the nozzle widens toward the exit) (m)
roof_r_at = lambda r: x_in_top - (r - r_in_rib_top) * roof_tan # roof height above radius r (m)
pass_reaches_manifold = (r_pass - d_pass / 2 >= r_in_rib_floor and r_pass + d_pass / 2 <= r_in_floor_outer
                         and roof_r_at(r_pass + d_pass / 2) >= x_in_floor + 1e-3) # passage pokes 0.7 mm in, roof above it
# Outlet port and passage
xs_op = np.linspace(x_out_top, x_op, 200)
metal_op_to_channels = np.min(r_op_pass - d_op_pass / 2 - ribtop_r(xs_op)) # passage to the rib tops (channels) under it (m)
op_pass_in_manifold = r_out_inner <= r_op_pass - d_op_pass / 2 and r_op_pass + d_op_pass / 2 <= r_out_outer
metal_op_port_to_pass = (r_op_top - d.get("orb4_tap_depth_m")) - (r_op_pass + d_op_pass / 2) # port bottom to passage top (m)
roof_x_at = lambda r: x_out_floor + (r_out_outer - r) * roof_tan # outlet roof height at radius r (m)
op_boss_clears_manifold = x_op - op_r >= roof_x_at(r_op_bottom) # boss bottom edge stays above the manifold roof
gap_op_to_bolts = r_bc * np.sin(np.radians(180 / n_bolts)) - d.get("bolt_head_d_m") / 2 - op_r # sideways, boss to the nearest bolt head (m)

print("-- Coolant Manifolds --")
print(f"Ethanol density in / out:      {rho_in:.0f} / {rho_out:.0f} kg/m^3")
print(f"Inlet manifold area:           {A_in_manifold*1e6:.1f} mm^2, velocity {v_in_manifold:.2f} m/s (channels there {v_ch_in:.1f} m/s, ratio {v_ch_in/v_in_manifold:.0f}x)")
print(f"Inlet port drill:              velocity {v_port:.2f} m/s")
print(f"Outlet manifold area:          {A_out_manifold*1e6:.1f} mm^2, velocity {v_out_manifold:.2f} m/s (channels there {v_ch_out:.1f} m/s, ratio {v_ch_out/v_out_manifold:.0f}x)")
print(f"Outlet port passage:           velocity {v_op_pass:.2f} m/s")
print(f"Pressure change round manifold: inlet {head_in*100:.2f} %, outlet {head_out*100:.2f} % of the jacket drop (even flow: < 1 %)")

print("\n-- Flange Bolts --")
print(f"Separation load:               {F_sep/1e3:.2f} kN (chamber gas out to the outer seal)")
print(f"Bolts:                         {n_bolts} x {d.get('bolt_size')}, preload {preload/1e3:.2f} kN each (75% of proof)")
print(f"Joint safety factor:           {joint_SF:.2f} (need >= {d.get('joint_safety_factor'):.1f})")
print(f"Flange bending stress:         {sigma_flange/1e6:.0f} MPa (allow {sigma_allow/1e6:.0f}), needs t >= {t_flange_needed*1e3:.1f} mm, have {t_flange*1e3:.1f} mm")

print("\n-- Geometry Checks --")
print(f"Metal, groove -> manifold:     {metal_groove_to_manifold*1e3:.2f} mm")
print(f"Metal, outlet roof -> outside: {metal_outlet_roof*1e3:.2f} mm (at the flange-top corner)")
print(f"Metal, inlet manifold -> exit: {metal_inlet_top*1e3:.2f} mm, -> collar outside {metal_inlet_outer*1e3:.2f} mm")
print(f"Channels open into manifolds:  {'yes' if inlet_opens_channels and rib_cut > 0 else 'NO'} (inner edge {rib_cut*1e3:.2f} mm into the rib tops)")
print(f"Flange covers jacket:          jacket OD at face {2*jacket_r(0)*1e3:.1f} mm < flange OD {d.get('flange_outer_d_m')*1e3:.0f} mm")
print(f"Port boss inside collar:       {'yes' if boss_fits else 'NO'} (boss x {(x_port-boss_r)*1e3:.2f} to {(x_port+boss_r)*1e3:.2f}, collar x {x_collar_start*1e3:.2f} to {L_total*1e3:.2f})")
print(f"Seal lands:                    bore -> inner groove {land_bore_inner*1e3:.2f} mm, inner -> outer groove {land_inner_outer*1e3:.2f} mm, outer groove -> bolt holes {land_outer_bolts*1e3:.2f} mm")
print(f"Metal under inlet passage:     {metal_under_passage*1e3:.2f} mm (jacket to passage, worst point from the port to the manifold)")
print(f"Passage opens into manifold:   {'yes' if pass_reaches_manifold else 'NO'}")
print(f"Outlet passage into manifold:  {'yes' if op_pass_in_manifold else 'NO'}, metal to channels {metal_op_to_channels*1e3:.2f} mm, port bottom to passage {metal_op_port_to_pass*1e3:.2f} mm")
print(f"Outlet boss clearances:        {'clear of' if op_boss_clears_manifold else 'FILLS'} the manifold roof, {gap_op_to_bolts*1e3:.1f} mm to the nearest bolt head")

problems = []
if joint_SF < d.get("joint_safety_factor"): problems.append("bolt preload too low")
if t_flange < t_flange_needed: problems.append("flange too thin")
if max(head_in, head_out) > 0.01: problems.append("manifold too small for even flow")
if min(metal_groove_to_manifold, metal_under_passage, metal_op_to_channels, metal_op_port_to_pass, land_bore_inner, land_inner_outer, land_outer_bolts) < 1e-3: problems.append("less than 1 mm metal between features")
if min(metal_outlet_roof, metal_inlet_top, metal_inlet_outer) < d.get("jacket_min_m"): problems.append("manifold wall thinner than jacket_min")
if not inlet_opens_channels: problems.append("inlet manifold misses some channel tops")
if not boss_fits: problems.append("port boss does not fit on the collar")
if not pass_reaches_manifold: problems.append("inlet passage misses the manifold")
if not op_pass_in_manifold: problems.append("outlet passage misses the manifold")
if not op_boss_clears_manifold: problems.append("outlet boss fills part of the manifold")
if gap_op_to_bolts < 2e-3: problems.append("outlet boss too close to the bolt heads")
print("\nProblems:                      " + (", ".join(problems) if problems else "none"))

# CAD dimensions (mm, x from injector face)
print("\n-- CAD Dimensions (mm) --")
print(f"Flange:          x 0 to {t_flange*1e3:.1f}, OD {d.get('flange_outer_d_m')*1e3:.1f}")
print(f"Bolt holes:      {n_bolts} x dia {d.get('bolt_clearance_d_m')*1e3:.1f} thru, on dia {d.get('bolt_circle_d_m')*1e3:.1f}, clock 22.5 + 45n")
print(f"Inner O-ring:    groove ID {d.get('inner_groove_id_m')*1e3:.2f}, OD {2*r_inner_seal*1e3:.2f}, depth {d.get('oring_groove_depth_m')*1e3:.2f} (on flange face)")
print(f"Outer O-ring:    groove ID {d.get('outer_groove_id_m')*1e3:.2f}, OD {2*r_outer_seal*1e3:.2f}, depth {d.get('oring_groove_depth_m')*1e3:.2f}")
fmt = lambda poly: ", ".join(f"({x*1e3:.2f}, {r*1e3:.2f})" for x, r in poly)
print(f"Outlet manifold: revolve (x, r): {fmt(out_poly)}")
print(f"Outlet port:     boss dia {2*op_r*1e3:.0f} at x {x_op*1e3:.2f} clock 0, from r {r_op_bottom*1e3:.2f} to flat top at r {r_op_top*1e3:.0f}; {d.get('outlet_port')}")
print(f"Outlet passage:  dia {d_op_pass*1e3:.0f} radial from port down to r {(r_op_pass - d_op_pass/2)*1e3:.1f} at x {x_op*1e3:.2f}, then axial at r {r_op_pass*1e3:.1f} from x {x_op*1e3:.2f} to {x_out_floor*1e3:.2f} (manifold floor)")
print(f"Channels:        start at x {x_ch_start*1e3:.1f} (under the outlet manifold), end at x {x_in_top*1e3:.2f} (top of the inlet manifold)")
print(f"Inlet collar:    x {x_collar_start*1e3:.2f} to {L_total*1e3:.2f}, OD {2*d.get('inlet_collar_outer_r_m')*1e3:.1f}")
print(f"Inlet manifold:  revolve (x, r): {fmt(in_poly)}")
print(f"Inlet port:      boss dia {d.get('inlet_port_boss_d_m')*1e3:.0f} at x {x_port*1e3:.2f}, flat top at r {r_boss_top*1e3:.0f}; {d.get('inlet_port')}")
print(f"Inlet passage:   dia {d_pass*1e3:.0f} radial from port down to r {r_pass*1e3:.1f} at x {x_port*1e3:.2f}, then axial at r {r_pass*1e3:.1f} from x {x_port*1e3:.2f} to {x_in_floor*1e3 + 0.7:.2f} (0.7 into the manifold floor)")
print(f"Gas side:        machining stock {d.get('machining_stock_m')*1e3:.1f}; nozzle exit at x {L_total*1e3:.2f}; chamber bore dia {2*R_chamber*1e3:.2f} (final)")

d.save_many({
    "inlet_manifold_area_mm2": (A_in_manifold*1e6, "mm^2"),
    "outlet_manifold_area_mm2": (A_out_manifold*1e6, "mm^2"),
    "x_outlet_manifold_top_mm": (x_out_top*1e3, "mm"),
    "x_inlet_manifold_floor_mm": (x_in_floor*1e3, "mm"),
    "x_channel_start_mm": (x_ch_start*1e3, "mm"),
    "x_channel_end_mm": (x_in_top*1e3, "mm"),
    "flange_separation_load_N": (F_sep, "N"),
    "flange_joint_SF": (joint_SF, "-"),
    "flange_bending_stress_Pa": (sigma_flange, "Pa"),
}, "chamber_hardware.py")
