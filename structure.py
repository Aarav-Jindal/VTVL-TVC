# Structure check: hand formulas for every part that carries load, at liftoff, full gimbal tilt with a gust, landing and
# handling. Prints a safety factor (SF) for each check against its minimum and saves the lowest. The checks, load cases and
# minimums were set on 2026-11-07. The closest parts get FEA separately.
# Run last, after plumbing_layout.py: it reads design_values.json, parts_list.csv and the layout stations.
# Checked by other scripts, not here: tank walls, caps and cap screws (feed_system.py), engine flanges and bolts
# (chamber_hardware.py, injector_hardware.py), gimbal ring, pins and arms (gimbal.py).
# Units: mm, N, MPa (N/mm^2), kg. u = mm above the injector face (the CAD's frame), up is +u.
import contextlib, csv, io, math, os, sys
import numpy as np
import design as d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1, os.path.join(HERE, "CAD", "scripts"))
with contextlib.redirect_stdout(io.StringIO()):
    import plumbing_routes as PR, mount_layout as ML, avionics_layout as AV, leg_layout as LL

g = 9.80665
SF_METAL, SF_COMP, SF_HW = d.get("structure_SF_metal_yield"), d.get("structure_SF_composite_ult"), d.get("structure_SF_screw_proof")

# ---- Materials, MPa (2026-11-07). [est] = typical for the material, not published for the exact part ----
# 6061-T6 [aerospacemetals.com ASM sheet]: bearing at e/D 2. Fbry at e/D 1.5 is the MIL-HDBK-5 value, not confirmed online [est]
AL = dict(Fty=276.0, Fsy=0.577*276.0, Fbry15=345.0, Fbry20=386.0, E=68900.0)
# G12: Wildman publishes no strengths. Typical filament-wound glass-epoxy tube [est], Gund F2016: 138 axial tension and compression,
# 48 shear. Bearing isn't published: taken as the compressive strength [est]. Axial modulus isn't published: 10 GPa [est, typical +/-55 deg wind]
G12 = dict(Fcu=138.0, Ftu=138.0, Fsu=48.0, Fbru=138.0, E=10000.0)
# G10 [Norplex NP510A]: tension 276, flexural 359 crosswise, punch shear 148, modulus 17.9 GPa. Bearing isn't published: taken as the
# edgewise compressive strength, 241 [Accurate Plastics minimum] [est]. Threads tapped edgewise in G10 have no published data: test a sample
G10 = dict(Ftu=276.0, Ffu=359.0, Fsu=148.0, Fbru=241.0, E=17900.0)
PETG_INSERT_N = 1167.0                                                       # Ruthex M3 insert in PETG, pull-out [CNC Kitchen]; M4 should be stronger [est]
SCREW = dict(proof=970.0, As=8.78, A_minor=7.75)                             # M4 class 12.9 [ISO 898-1]: grade not picked yet (placeholder)
M6_SHANK = math.pi/4*6.0**2                                                  # hinge pin: an M6x60 screw, plain shank in the lugs
CABLE_BREAK_N = 4090.0                                                       # 3/32 in 7x19 302/304, 920 lbf minimum [Loos, Lawson]
SS316 = dict(Fty=205.0, E=193000.0)                                          # ASTM A269 TP316 annealed, spec minimum
CARBON_E = (70000.0, 90000.0)                                                # roll-wrapped strut, 70 [est] to 90 [listed for the 50.3 mm tube]
CARBON_FLEX = 600.0                                                          # strut tube bending strength [est], only used for a reported bound

# ---- Geometry, from the layout data ----
R_O, R_I = PR.R_WALL_O, PR.R_WALL_I; T_W = R_O - R_I                         # body tube, 2.4 mm wall
D4 = 4.0                                                                     # M4 shank
L_GAS = d.get("L_gas_total_mm")
u_of = lambda x_m: x_m*1e3 - L_GAS                                           # m above the nozzle exit -> u
U_BOT, U_TOP = AV.U_BOT, AV.U_TOP
U_PIVOT = d.get("gimbal_pivot_above_face_m")*1e3
U_MOUNT, U_RING_LEG = LL.U_MOUNT, LL.U_RING
CLOCKS = LL.CLOCKS

# ---- Masses along the body: laid out like vehicle_layout.py, then dry parts scaled up to the budget's margins ----
gap = d.get("section_gap_m"); cap = d.get("bulkhead_length_m")
order = [("engine_bay", d.get("engine_bay_len_m")), ("n2o_tank", d.get("L_tank_N2O_m") + cap), ("eth_tank", d.get("L_tank_ethanol_m") + cap),
         ("he_valves", d.get("He_valve_stack_len_m")), ("he_bottle", d.get("He_bottle_len_m")), ("avionics", d.get("avionics_bay_len_m"))]
sec = {}; x = 0.0
for i, (name, length) in enumerate(order):
    if i and name != "he_bottle": x += d.get("n2o_eth_gap_m") if name == "eth_tank" else gap
    sec[name] = (x, x + length); x += length
mid = lambda s: 0.5*(sec[s][0] + sec[s][1])
x_face = L_GAS*1e-3
spots = {"avionics": mid("avionics"), "he_valves": mid("he_valves"), "he_bottle": mid("he_bottle"),
         "eth_top": sec["eth_tank"][1] - cap/4, "eth_bottom": sec["eth_tank"][0] + cap/4,
         "n2o_top": sec["n2o_tank"][1] - cap/4, "n2o_bottom": sec["n2o_tank"][0] + cap/4,
         "engine_bay": 0.5*(x_face + d.get("inj_hub_thickness_m") + d.get("pintle_head_len_m") + sec["engine_bay"][1]),
         "engine": x_face - d.get("x_engine_cg_mm")*1e-3, "engine_top": x_face + d.get("inj_hub_thickness_m") + 0.02,
         "gimbal": x_face + d.get("gimbal_pivot_above_face_m"), "legs": -d.get("leg_cg_below_exit_m")}
with open(os.path.join(HERE, "parts_list.csv"), newline="", encoding="utf-8") as f:
    rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]
kg = lambda r: int(r["qty"])*float(r["each_g"])/1000
group = lambda w: sum(kg(r) for r in rows if r["where"] == w)
dry = []; spread = d.get("m_body_tube_kg")                                   # (name, kg, u); spread = kg along the whole tube
for r in rows:
    w = r["where"]
    if w == "body" and r["item"].startswith("Body tube joints"): dry += [(r["item"], kg(r)/2, u) for u in AV.U_JOINT]
    elif w == "body": spread += kg(r)
    elif w == "rings": dry += [(r["item"], kg(r)/4, u) for u in ML.U_TANK_RINGS]
    else: dry.append((r["item"], kg(r), u_of(spots[w])))
dry += [("Engine (CAD)", d.get("m_engine_cad_kg"), u_of(spots["engine"])), ("Gimbal (CAD)", d.get("m_gimbal_cad_kg"), U_PIVOT),
        ("Nitrous tank", d.get("m_tank_N2O_kg"), 0.5*(PR.u_nb + PR.u_nt)), ("Ethanol tank", d.get("m_tank_ethanol_kg"), 0.5*(PR.u_eb + PR.u_et))]
m_bottle = d.get("He_bottle_mass_kg")/2
U_BOTTLES = (PR.u_he + ML.LEN/2, PR.u_he2 + ML.LEN/2)
dry += [("Helium bottle", m_bottle, u) for u in U_BOTTLES]
N_SEG = 44; seg = np.linspace(U_BOT, U_TOP, N_SEG + 1)
dry += [("Body tube and spread parts", spread/N_SEG, 0.5*(a + b)) for a, b in zip(seg, seg[1:])]
m_dry_list = sum(m for _, m, _ in dry)
K_MARGIN = d.get("m_dry_kg")/m_dry_list                                     # every dry part gets the budget's average margin
dry = [(n, m*K_MARGIN, u) for n, m, u in dry]
r_in = d.get("tank_outer_radius_m") - d.get("tank_wall_m"); A_in = math.pi*r_in**2
def liquid(name, tank, m, rho):
    h = m/(rho*A_in); return (name, m, u_of(sec[tank][0] + cap/2 + h/2))
m_He = d.get("m_He_kg")
full = dry + [liquid("Nitrous", "n2o_tank", d.get("m_N2O_kg"), d.get("rho_N2O")), liquid("Ethanol", "eth_tank", d.get("m_ethanol_kg"), d.get("rho_ethanol"))] \
       + [("Helium", m_He/2, u) for u in U_BOTTLES]
empty = dry + [("Helium", m_He/2, 0.5*(PR.u_nb + PR.u_nt)), ("Helium", m_He/2, 0.5*(PR.u_eb + PR.u_et))]
m_legs = sum(m for n, m, _ in dry if n.startswith("Landing legs"))
empty_body = [l for l in empty if not l[0].startswith("Landing legs")]     # at touchdown the legs stop themselves; the body is what the struts stop
M_FULL, M_EMPTY = sum(l[1] for l in full), sum(l[1] for l in empty)

def internal(lumps, loads, cuts):
    """Free body: lumps (name, kg, u) on the axis, external non-gravity loads (u, Fa up, Fl sideways, Mp N mm).
    Returns the felt axial acceleration (m/s^2) and, at each cut, what the part above passes to the part below:
    N (+ = tension), V and M (N mm)."""
    m = np.array([l[1] for l in lumps]); u = np.array([l[2] for l in lumps]); M = m.sum(); uc = (m*u).sum()/M
    I = (m*(u - uc)**2).sum()
    aa = sum(l[1] for l in loads)/M
    alpha = sum((l[0] - uc)*l[2] + l[3] for l in loads)/I
    al = sum(l[2] for l in loads)/M + alpha*(u - uc)
    out = []
    for c in cuts:
        b = u < c; ext = [l for l in loads if l[0] < c]
        N = (m[b]*aa).sum() - sum(l[1] for l in ext)
        V = (m[b]*al[b]).sum() - sum(l[2] for l in ext)
        Mc = (m[b]*al[b]*(u[b] - c)).sum() - sum((l[0] - c)*l[2] + l[3] for l in ext)
        out.append((N, V, Mc))
    return aa, out

# ---- Leg forces at full crush: the cartridge caps the strut at leg_crush_force_N; the cables and the ground balance it at the foot ----
S = d.get("leg_crush_force_N")
def leg_forces(c):
    node = LL.P(LL.U_FOOT, LL.R_FOOT, c); pin = LL.strut_pin(c)
    e_s = (node - pin)/np.linalg.norm(node - pin)                            # strut pushes the foot along this
    e_c = [(LL.anchor_pin(c, s) - node)/np.linalg.norm(LL.anchor_pin(c, s) - node) for s in (-1, 1)]
    up = np.array([-1.0, 0, 0])
    A = np.column_stack([e_c[0] + e_c[1], up]); Cc, R = np.linalg.lstsq(A, -S*e_s, rcond=None)[0]
    return dict(R=R, C=Cc, pin_force=-S*e_s, pin=pin, anchor_force={s: Cc*LL.cable_dir(c, s) for s in (-1, 1)},
                anchor={s: LL.anchor_pin(c, s) for s in (-1, 1)})
LEG = leg_forces(CLOCKS[0])
ax = lambda v: -v[0]                                                         # axial (up) part of a rocket-frame vector
radial = lambda v, c: float(np.dot(v, LL.rad(c)))
R_FOOT_N, C_CABLE = LEG["R"], LEG["C"]
N_LAND = 4*R_FOOT_N/(sum(l[1] for l in empty_body)*g)                       # felt g on the body, all 4 legs crushing
S_UP = ax(LEG["pin_force"]); S_IN = -radial(LEG["pin_force"], CLOCKS[0])
MU = 0.5                                                                     # foot friction, as in the landing simulation
N_LIFT = d.get("thrust_N")/(M_FULL*g)

# ---- Load cases along the body ----
T = d.get("thrust_N"); TILT = math.radians(d.get("gimbal_range_deg"))
gust = 0.5*d.get("air_density_kgm3")*d.get("gust_mps")**2*d.get("Cd_side")*d.get("body_od_m")*d.get("L_body_tube_m")
U_GUST = u_of(d.get("x_body_mid_m"))
def leg_loads(c, one=False):
    """External loads on the body from one leg in the plane of the page (its clock at +y), or all 4 (symmetric: axial only)."""
    lf = leg_forces(c); out = []
    for pt, F in [(lf["pin"], lf["pin_force"])] + [(lf["anchor"][s], lf["anchor_force"][s]) for s in (-1, 1)]:
        y = radial(pt, c); Fa = ax(F)
        out.append((-pt[0], 4*Fa if not one else Fa, 0.0 if not one else radial(F, c), 0.0 if not one else -y*Fa))
    return out
cases = {"liftoff": (full, [(U_PIVOT, T, 0.0, 0.0)])}
for name, lumps in (("full", full), ("empty", empty)):
    for sgn in (1, -1):
        cases[f"tilt {name}, gust {'+' if sgn > 0 else '-'}"] = (lumps, [(U_PIVOT, T*math.cos(TILT), T*math.sin(TILT), 0.0), (U_GUST, 0.0, sgn*gust, 0.0)])
cases["landing, 4 legs"] = (empty_body, leg_loads(CLOCKS[0]))
cases["landing, 1 leg first"] = (empty_body, leg_loads(CLOCKS[0], one=True))
U_CRADLE = (U_MOUNT, AV.U_BAY_RING)                                         # lying on its side: cradles under the engine mount and the bay ring
W_DRY = sum(l[1] for l in empty)*g; u_cg_dry = sum(l[1]*l[2] for l in empty)/sum(l[1] for l in empty)
R2 = W_DRY*(u_cg_dry - U_CRADLE[0])/(U_CRADLE[1] - U_CRADLE[0]); R1 = W_DRY - R2
cases["handling"] = (empty, [(U_CRADLE[0], 0.0, R1, 0.0), (U_CRADLE[1], 0.0, R2, 0.0)])
CUTS = np.linspace(U_BOT + 1, U_TOP - 1, 400)
loads_at = {k: internal(v[0], v[1], list(CUTS) + list(AV.U_JOINT)) for k, v in cases.items()}

# ---- Checks ----
results = []
def check(part, what, case, sf, sf_min, basis):
    results.append((part, what, case, sf, sf_min, basis))
def joint_screw_force(N, V, M, n=8, r=R_I):
    return abs(N)/n + 2*abs(M)/(n*r) + 2*abs(V)/n                           # worst screw of a ring of n: axial share + moment + shear
tau_screw = 0.577*SCREW["proof"]                                             # shear allowable from the proof stress (von Mises)

# Body tube, everywhere: axial + bending against compression (or shell buckling) and tension
A_t = math.pi*(R_O**2 - R_I**2); Z_t = math.pi/32*((2*R_O)**4 - (2*R_I)**4)/(2*R_O)
R_m = 0.5*(R_O + R_I)
knock = 1 - 0.901*(1 - math.exp(-math.sqrt(R_m/T_W)/16))                    # NASA SP-8007 knockdown
sig_cr = knock*0.605*G12["E"]*T_W/R_m if G12["E"] else 0.0                 # isotropic formula with the axial modulus: rough for a wound tube
worst_c = worst_t = (np.inf, "")
for k, (aa, out) in loads_at.items():
    for (N, V, M), u in zip(out[:len(CUTS)], CUTS):
        sc = -N/A_t + abs(M)/Z_t; st = N/A_t + abs(M)/Z_t
        if sc > 0: worst_c = min(worst_c, (min(G12["Fcu"], sig_cr)/sc, f"{k}, u {u:.0f}"))
        if st > 0: worst_t = min(worst_t, (G12["Ftu"]/st, f"{k}, u {u:.0f}"))
check("Body tube", "compression + bending (or buckling)", worst_c[1], worst_c[0], SF_COMP, "G12 ult")
check("Body tube", "tension + bending", worst_t[1], worst_t[0], SF_COMP, "G12 ult")

# Joint rings: 8 M4 each side, through the tube into the 3 mm ring wall
e_joint = AV.ROW - AV.GAP/2; t_jring = (2*AV.R_FIT - AV.RING_ID)/2
for j, uj in enumerate(AV.U_JOINT):
    Fw, kw = 0.0, ""
    for k, (aa, out) in loads_at.items():
        F = joint_screw_force(*out[len(CUTS) + j])
        if F > Fw: Fw, kw = F, k
    part = f"Joint {j + 1} (u {uj:.0f})"
    check(part, f"M4 shear ({Fw:.0f} N, worst screw)", kw, tau_screw/(Fw/SCREW["A_minor"]), SF_HW, "12.9 proof")
    check(part, "tube bearing", kw, G12["Fbru"]/(Fw/(D4*T_W)), SF_COMP, "G12 ult")
    check(part, f"tube tear-out to its end ({e_joint:.1f} mm)", kw, G12["Fsu"]/(Fw/(2*(e_joint - D4/2)*T_W)), SF_COMP, "G12 ult")
    check(part, f"ring bearing ({t_jring:.1f} mm wall)", kw, AL["Fbry20"]/(Fw/(D4*t_jring)), SF_METAL, "6061 yield")

# Engine mount: 8 M4 into the gimbal mount. Thrust at full tilt; at landing the cable lugs sit on these screws
e_mount = U_MOUNT - U_BOT
F_lift = T*math.cos(TILT)/8 + 2*T*math.sin(TILT)/8
m_swing = sum(m for n, m, u in dry if u < U_PIVOT + 30)                     # engine, gimbal and what's on them
F_eng_land = m_swing*N_LAND*g/8
cab = LEG["anchor_force"][1]; c_an = CLOCKS[0] + LL.ANCHOR_OFFSET
F_tens = max(radial(cab, c_an), 0.0); F_sh = math.hypot(ax(cab), float(np.dot(cab, LL.tan(c_an)))) + F_eng_land
sf_comb = 1/math.hypot((F_tens/SCREW["As"])/SCREW["proof"], (F_sh/SCREW["A_minor"])/tau_screw)
check("Engine mount screws", f"M4 shear, thrust at full tilt ({F_lift:.0f} N)", "tilt", tau_screw/(F_lift/SCREW["A_minor"]), SF_HW, "12.9 proof")
check("Engine mount screws", f"M4 tension + shear under a cable lug ({F_tens:.0f} + {F_sh:.0f} N)", "landing", sf_comb, SF_HW, "12.9 proof")
L_e_mount = 12.0 - 3.0 - T_W                                                 # M4x12 through the 3 mm lug floor and the tube
check("Engine mount screws", f"thread pull-out of the mount ({L_e_mount:.1f} mm engaged)", "landing", AL["Fsy"]/(F_tens/(0.5*math.pi*D4*L_e_mount)), SF_METAL, "6061 yield")
Fm = max(F_lift, F_sh)
check("Engine mount screws", "tube bearing", "tilt" if F_lift > F_sh else "landing", G12["Fbru"]/(Fm/(D4*T_W)), SF_COMP, "G12 ult")
check("Engine mount screws", f"tube tear-out to its end ({e_mount:.0f} mm)", "tilt" if F_lift > F_sh else "landing", G12["Fsu"]/(Fm/(2*(e_mount - D4/2)*T_W)), SF_COMP, "G12 ult")

# Rings that carry a mass: tank bottoms, bottle seats, the bay ring. Load = what sits on the ring x the worse of liftoff and landing.
# Parts get the budget's margins: 1.1 on known masses, 1.25 on estimates
k_est, k_known = 1.25, 1.1
def ring_load(m_full, m_empty):
    lift, land = m_full*N_LIFT*g, m_empty*N_LAND*g
    return (lift, "liftoff") if lift > land else (land, "landing")
n2o_hw = group("n2o_top") + group("n2o_bottom"); eth_hw = group("eth_top") + group("eth_bottom")
reg_one = kg(next(r for r in rows if r["item"].startswith("Regulators")))/2
fit = kg(next(r for r in rows if r["item"].startswith("Helium end fittings")))
stack_up = reg_one + fit/2                                                   # upper bottle: its regulator and half the fittings [est split]
stack_low = group("he_valves") - stack_up                                    # lower bottle: the rest of the stack hangs off its neck
m_t_n, m_t_e = d.get("m_tank_N2O_kg")*k_known, d.get("m_tank_ethanol_kg")*k_known
rings = [
    ("Nitrous tank ring (bottom)", "al", ring_load(m_t_n + n2o_hw*k_est + d.get("m_N2O_kg"), m_t_n + n2o_hw*k_est + m_He/2), 79.4, ML.R_RING, 76.35, 17.0),
    ("Ethanol tank ring (bottom)", "al", ring_load(m_t_e + eth_hw*k_est + d.get("m_ethanol_kg"), m_t_e + eth_hw*k_est + m_He/2), 79.4, ML.R_RING, 76.35, 17.0),
    ("Lower bottle seat ring", "g10", ring_load(m_bottle*k_known + stack_low*k_est + m_He/2, m_bottle*k_known + stack_low*k_est), ML.R_BORE_MID + 1, ML.R_RING, ML.R_BORE_MID, ML.ND),
    ("Upper bottle seat ring", "g10", ring_load(m_bottle*k_known + stack_up*k_est + m_He/2, m_bottle*k_known + stack_up*k_est), ML.R_BORE_MID + 1, ML.R_RING, ML.R_BORE_MID, ML.ND),
    ("Bay ring", "g10", ring_load(0, (group("avionics") - 0.237)*k_est), AV.TRAY_SCREW_R, AV.R_FIT, AV.BAY_RING_ID/2, 0.0)]
for part, mat, (F, case), r_load, r_out, r_in_, notch in rings:
    t_r = 6.0 if mat == "al" else ML.T_RING; Fs = F/8; L_e = 8.0 - T_W        # M4x8 through the tube: 5.6 mm in the ring
    check(part, f"M4 shear ({Fs:.0f} N each)", case, tau_screw/(Fs/SCREW["A_minor"]), SF_HW, "12.9 proof")
    check(part, "tube bearing", case, G12["Fbru"]/(Fs/(D4*T_W)), SF_COMP, "G12 ult")
    if mat == "al": check(part, "ring bearing on the thread", case, AL["Fbry20"]/(Fs/(D4*L_e)), SF_METAL, "6061 yield")
    else: check(part, "G10 bearing on the thread", case, G10["Fbru"]/(Fs/(D4*L_e)), SF_COMP, "G10 ult")
    w = F/(2*math.pi*r_load); L_span = 2*math.pi*0.5*(r_out + r_in_)/8     # load spread round the ring, carried to 8 screws
    b_net = (r_out - r_in_) - notch
    sig = 6*(w*L_span**2/12)/(b_net*t_r**2) + 6*(w*(r_out - r_load))/(t_r**2) # span between screws (fixed ends) + the ring as a radial cantilever
    allow, sfm, basis = (AL["Fty"], SF_METAL, "6061 yield") if mat == "al" else (G10["Ffu"], SF_COMP, "G10 ult")
    check(part, f"ring bending between screws ({b_net:.1f} mm {'at a notch' if notch else 'section'})", case, allow/sig, sfm, basis)

# Legs: strut, hinge pin and lugs, the bracket's 2 screws into the leg ring, the leg ring itself, the cables
I_str = math.pi/64*(LL.R_TUBE_O*2)**4 - math.pi/64*(LL.R_TUBE_I*2)**4
for E in CARBON_E:
    check("Leg strut", f"buckling at the crush force (E {E/1e3:.0f} GPa)", "landing", math.pi**2*E*I_str/LL.L_STRUT**2/S, SF_COMP, "Euler, pinned")
check("Leg hinge pin (M6 shank)", "double shear", "landing", tau_screw/(S/(2*M6_SHANK)), SF_HW, "12.9 proof (grade not picked)")
check("Leg hinge bracket", "lug bearing (2 x 5 mm)", "landing", AL["Fbry15"]/(S/(2*6.0*5.0)), SF_METAL, "6061 yield")
check("Leg hinge bracket", "lug tear-out (11 mm edge)", "landing", AL["Fsy"]/(S/(2*2*(11.0 - 3.1)*5.0)), SF_METAL, "6061 yield")
check("Leg bracket screws (2 per leg)", f"M4 shear ({S_UP/2:.0f} N each)", "landing", tau_screw/(S_UP/2/SCREW["A_minor"]), SF_HW, "12.9 proof")
check("Leg bracket screws (2 per leg)", "tube bearing", "landing", G12["Fbru"]/(S_UP/2/(D4*T_W)), SF_COMP, "G12 ult")
check("Leg bracket screws (2 per leg)", "leg ring bearing on the thread", "landing", AL["Fbry20"]/(S_UP/2/(D4*(8.0 - T_W))), SF_METAL, "6061 yield")
th = math.pi/4; R_lr = 0.5*(ML.R_RING + 76.35); t_lr = 6.0; b_lr = (ML.R_RING - 76.35) - 17.0 # Tank_Ring reused: 6 mm, notches 17 deep
def ring_moment(phi):
    """Ring with 4 equal inward point loads (Roark): bending moment and hoop compression phi rad from a load."""
    a = th - phi                                                             # from the midpoint between loads
    return S_IN*R_lr/2*(1/th - math.cos(a)/math.sin(th)), S_IN*math.cos(a)/(2*math.sin(th))
NOTCH_CLOCKS = (30.0, 60.0, 210.0, 240.0)                                    # Tank_Ring's notches; the leg ring is mated square to the body's planes.
phi_n = min(min(abs((nc - lc + 180) % 360 - 180) for lc in CLOCKS) for nc in NOTCH_CLOCKS) # turned by any 45 deg step, the nearest notch is still 7.5 deg from a leg
M_lr, N_lr = ring_moment(math.radians(phi_n))
sig_lr = abs(M_lr)/(t_lr*b_lr**2/6) + N_lr/(t_lr*b_lr)
check("Leg ring", f"bending at the notch {phi_n:.1f} deg from a leg ({S_IN:.0f} N in at each leg)", "landing", AL["Fty"]/sig_lr, SF_METAL, "6061 yield, point loads")
check("Leg cables", f"tension ({C_CABLE:.0f} N each)", "landing", CABLE_BREAK_N/C_CABLE, SF_HW, "min breaking (screw rule)")

# Foot caught sideways (reported only): a friction-limited side load at the foot. How it splits between the strut (bending off the
# clevis) and the cable V depends on stiffnesses hand formulas can't pin down, so the two bounds are printed: all on the strut, all on the cable
F_side = MU*R_FOOT_N
M_root = F_side*LL.L_STRUT
SF_side_strut = CARBON_FLEX/(M_root/(I_str/LL.R_TUBE_O)) if CARBON_FLEX else 0.0
SF_side_lug = AL["Fbry15"]/(M_root/(5.0 + 40.0)/(6.0*5.0)) if AL["Fbry15"] else 0.0
beta = abs(float(np.dot(LL.cable_dir(CLOCKS[0], 1), LL.tan(CLOCKS[0]))))  # how much of the tight cable points sideways
T_side_cable = C_CABLE + F_side/beta

# Hanging parts: relief valves on 1/4 in tube, landing (down) plus friction-limited side load
I_14 = math.pi/64*(6.35**4 - (6.35 - 2*0.889)**4); Z_14 = I_14/(6.35/2)
m_rv = 0.400*k_est; W_rv = m_rv*g
rv = [("RV-2 on the nitrous top cap", PR.R_RV - PR.R_PORT, PR.u_rv2 - PR.u_ct_n + PR.L_RV/2),
      ("RV-3 on the ethanol top cap", PR.R_RV - PR.R_PORT, PR.u_rv3 - PR.u_ct_e + PR.L_RV/2),
      ("RV-1 on the helium stack", PR.R_RV - (PR.PORT_BOSS + PR.H_T4 - PR.D14), PR.u_rv1 - PR.u_side_rv + PR.L_RV/2)]
for part, arm, h in rv:
    sig = W_rv*(N_LAND*arm + MU*N_LAND*h)/Z_14
    check(part, f"tube bending at its root ({arm:.0f} mm out, {h:.0f} mm up)", "landing", SS316["Fty"]/sig, SF_METAL, "316 annealed yield")

# Fasteners at the top: dome on 8 M4 into heat-set inserts (side load at landing)
F_dome = 0.237*k_known*MU*N_LAND*g/8
check("Dome inserts (8)", f"insert load ({F_dome:.1f} N each)", "landing", PETG_INSERT_N/F_dome if PETG_INSERT_N else 0.0, SF_COMP, "insert pull-out [est]")

# ---- Reported only: helium stack moment on the bottle neck, bending mode, hanging parts' frequencies ----
m_stack = stack_low*k_est; h_stack = 110.0                                  # CG below the neck [est]
M_neck = m_stack*MU*N_LAND*g*h_stack

def bending_modes(lumps, with_tanks):
    nodes = np.linspace(U_BOT, U_TOP, 89); n = len(nodes); K = np.zeros((2*n, 2*n)); Mm = np.zeros((2*n, 2*n))
    EI_t = G12["E"]*math.pi/64*((2*R_O)**4 - (2*R_I)**4)
    r_t = d.get("tank_outer_radius_m")*1e3; I_tank = math.pi/64*((2*r_t)**4 - (2*(r_t - d.get("tank_wall_m")*1e3))**4)
    spans = [(ML.U_TANK_RINGS[0], ML.U_TANK_RINGS[1]), (ML.U_TANK_RINGS[2], ML.U_TANK_RINGS[3])]
    for i in range(n - 1):
        L = nodes[i + 1] - nodes[i]; um = 0.5*(nodes[i] + nodes[i + 1])
        EI = EI_t + (AL["E"]*I_tank if with_tanks and any(a <= um <= b for a, b in spans) else 0.0)
        k = EI/L**3*np.array([[12, 6*L, -12, 6*L], [6*L, 4*L*L, -6*L, 2*L*L], [-12, -6*L, 12, -6*L], [6*L, 2*L*L, -6*L, 4*L*L]])
        idx = [2*i, 2*i + 1, 2*i + 2, 2*i + 3]; K[np.ix_(idx, idx)] += k
    for _, m, u in lumps:                                                   # each lump shared between its two nearest nodes, kg -> N s^2/mm
        i = int(np.clip(np.searchsorted(nodes, u) - 1, 0, n - 2)); f = np.clip((u - nodes[i])/(nodes[i + 1] - nodes[i]), 0, 1)
        Mm[2*i, 2*i] += (1 - f)*m*1e-3; Mm[2*i + 2, 2*i + 2] += f*m*1e-3
    md = np.diag(Mm).copy(); md[md < 1e-12] = 0.0
    mst = np.where(md > 0)[0]; slv = np.where(md == 0)[0]                  # lumped masses: condense out the massless DOFs (exact for lumps)
    K_red = K[np.ix_(mst, mst)] - K[np.ix_(mst, slv)] @ np.linalg.solve(K[np.ix_(slv, slv)], K[np.ix_(slv, mst)])
    Dh = np.diag(1/np.sqrt(md[mst]))
    w2 = np.sort(np.linalg.eigvalsh(Dh @ (0.5*(K_red + K_red.T)) @ Dh))
    assert abs(w2[1]) < 1e-3*w2[2], "the two rigid-body modes should come out near 0"
    return math.sqrt(w2[2])/(2*math.pi)                                     # first flexible mode, free-free

if __name__ == "__main__":
    print("-- Loads --")
    print(f"Masses: {M_FULL:.2f} kg full, {M_EMPTY:.2f} kg empty (budget {d.get('m_wet_kg'):.2f} and {d.get('m_dry_kg') + m_He:.2f}); dry parts x {K_MARGIN:.3f} for margin")
    print(f"Liftoff: {N_LIFT:.2f} g. Full tilt: {d.get('gimbal_range_deg'):.0f} deg, side force {T*math.sin(TILT):.0f} N; gust {gust:.1f} N at u {U_GUST:.0f}")
    print(f"Landing, every strut at {S:.0f} N: foot {R_FOOT_N:.0f} N up, each cable {C_CABLE:.0f} N, strut on the bracket {S_UP:.0f} N up and {S_IN:.0f} N in")
    print(f"  {N_LAND:.2f} g on the body (legs excluded); friction {MU} at the feet gives up to {MU*N_LAND:.2f} g sideways")
    print(f"Handling: cradles at u {U_CRADLE[0]:.0f} and {U_CRADLE[1]:.0f}, {R1:.0f} and {R2:.0f} N")
    for k, (aa, out) in loads_at.items():
        Ms = [abs(o[2]) for o in out[:len(CUTS)]]; i = int(np.argmax(Ms))
        Nc = min(o[0] for o in out[:len(CUTS)])
        print(f"  {k:26s} largest bending {Ms[i]/1e3:6.0f} N m at u {CUTS[i]:5.0f}; most compression {-Nc:6.0f} N")

    print(f"\n-- Safety factors (minimum: metal {SF_METAL:.0f} on yield, G12/G10/PETG {SF_COMP:.0f} on ultimate, screws and cables {SF_HW:.0f} on proof or breaking) --")
    last = None
    for part, what, case, sf, sfm, basis in results:
        if part != last: print(part); last = part
        flag = "  LOW" if sf < sfm else ""
        print(f"  {what:58s} {case:24s} SF {sf:7.2f}  (min {sfm:.0f}, {basis}){flag}")

    f_tube = bending_modes(full, False); f_tanks = bending_modes(full, True)
    print("\n-- Reported, no pass or fail --")
    print(f"First bending mode, full: {f_tube:.0f} Hz with the tube alone, {f_tanks:.0f} Hz counting the tanks as stiffening it (rough beam)")
    print(f"Helium stack on the lower bottle's neck: {m_stack:.2f} kg, {M_neck/1e3:.1f} N m at landing sideways (no published neck or regulator strength to compare)")
    print(f"Foot caught sideways, {F_side:.0f} N (friction {MU} x the foot load): if the strut takes it all, strut SF {SF_side_strut:.2f} and hinge lugs SF {SF_side_lug:.2f};"
          f" if the cables take it all, the tight one sees {T_side_cable:.0f} N (SF {CABLE_BREAK_N/T_side_cable:.2f} on breaking). The real split needs FEA or a bench test")
    for part, arm, h in rv:
        f_rv = math.sqrt(3*SS316["E"]*I_14/(h**3)/(m_rv*1e-3))/(2*math.pi)
        print(f"{part}: about {f_rv:.0f} Hz sideways on its tube (upper bound: ignores the bend's give)")

    low = [r for r in results if r[3] < r[4]]
    worst = min(results, key=lambda r: r[3]/r[4])
    print("\nProblems:                      " + ("none" if not low else "; ".join(f"{r[0]}: {r[1]} SF {r[3]:.2f} < {r[4]:.0f}" for r in low)))
    d.save_many({"structure_landing_g": (N_LAND, "g"),
                 "structure_min_SF_ratio": (worst[3]/worst[4], "-"),
                 "structure_min_SF_item": (f"{worst[0]}: {worst[1]}", "-"),
                 "structure_bending_Hz_tube": (f_tube, "Hz"),
                 "structure_bending_Hz_with_tanks": (f_tanks, "Hz")}, "structure.py")
