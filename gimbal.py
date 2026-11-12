# Gimbal (TVC) sizing: how far the engine must swing to hold the vehicle (misalignment, CG offset, a side gust, maneuvers),
# the mass and inertia that swing with it, the torque and speed the gimbal axis needs, and the servo linkage that gives both.
# Heights are measured up from the nozzle exit (x 0), as in vehicle_layout.py. Run after vehicle_layout.py.
import csv
import numpy as np
import design as d

g = 9.80665
mm = 1e-3
rad = np.radians
F_max = d.get("thrust_N")
rng = d.get("gimbal_range_deg"); use = d.get("gimbal_use_max")

# ---- What swings: the engine (CAD) + the parts mounted on it (parts_list.csv: engine, engine_top) ----
x_face = d.get("L_gas_total_mm")*mm
x_pivot = d.get("x_gimbal_pivot_m"); p = x_pivot - x_face # pivot height above the injector face
m_eng = d.get("m_engine_cad_kg"); x_eng = x_face - d.get("x_engine_cg_mm")*mm
I_piv = d.get("I_engine_face_kgm2") + m_eng*(x_eng - x_pivot)**2 - m_eng*(x_eng - x_face)**2 # parallel axis: face centre -> pivot
spot = {"engine": x_eng, "engine_top": x_face + d.get("inj_hub_thickness_m") + 0.02} # same spots as vehicle_layout.py
m_move, mx = m_eng, m_eng*x_eng
with open("parts_list.csv", newline="", encoding="utf-8") as f:
    for r in csv.DictReader(line for line in f if not line.startswith("#")):
        if r["where"] in spot:
            m = int(r["qty"])*float(r["each_g"])/1000; x = spot[r["where"]]
            m_move += m; mx += m*x; I_piv += m*(x - x_pivot)**2
mm_arm = d.get("m_gimbal_arm_kg"); r_arm = 0.5*(d.get("gimbal_arm_riser_r_in_m") + d.get("gimbal_arm_riser_r_out_m"))
I_piv += 3*mm_arm*(r_arm**2/2 + d.get("gimbal_arm_top_above_pivot_m")**2) # 3 arms swing with the engine, spread around it
m_move += 3*mm_arm; mx += 3*mm_arm*x_pivot
x_move = mx/m_move # center of mass of everything that swings
R_ring = 0.25*(d.get("gimbal_ring_id_m") + d.get("gimbal_ring_od_m"))
I_piv += d.get("m_gimbal_ring_kg")*R_ring**2/2 # the ring also swings about the outer axis (worst axis)

# ---- Angle needed to hold the vehicle at hover (thrust = weight), liftoff and landing ----
gust_F = 0.5*d.get("air_density_kgm3")*d.get("gust_mps")**2*d.get("Cd_side")*d.get("body_od_m")*d.get("L_body_tube_m")
x_cp = d.get("x_body_mid_m") # side force acts at the middle of the body tube
M_full = d.get("m_wet_kg"); M_empty = M_full - d.get("m_prop_kg") # mass loop values (with margin)
cases = {"liftoff": (M_full, d.get("x_cg_liftoff_m"), d.get("I_pitch_liftoff_kgm2")),
         "landing": (M_empty, d.get("x_cg_landing_m"), d.get("I_pitch_landing_kgm2"))}
angle = lambda torque, T, arm: np.degrees(np.arcsin(min(torque/(T*arm), 1.0))) # gimbal angle that gives this torque
need = {}
for case, (M, x_cg, I_veh) in cases.items():
    T = M*g; arm = x_cg - x_pivot
    parts = {"misalignment": d.get("thrust_misalign_deg"),
             "CG offset": angle(M*g*d.get("cg_offset_m"), T, arm),
             "gust": angle(gust_F*abs(x_cp - x_cg), T, arm),
             "maneuver": angle(I_veh*rad(d.get("maneuver_accel_dps2")), T, arm)}
    need[case] = (parts, sum(parts.values()), arm, T)
worst = max(need, key=lambda c: need[c][1])
angle_req = need[worst][1]

# ---- Torque and speed at the gimbal axis ----
A = rad(use*rng); w_bw = 2*np.pi*d.get("actuator_bandwidth_hz") # follow a sine of the usable amplitude at the bandwidth
torque = {"swinging inertia": I_piv*A*w_bw**2,
          "hose stiffness": d.get("hose_stiffness_Nm_per_deg")*rng,
          "engine weight off the pivot": m_move*g*abs(x_move - x_pivot)*np.sin(rad(rng)),
          "thrust line off the pivot": F_max*d.get("pivot_offset_m"),
          "friction": d.get("gimbal_friction_Nm")}
tq_req = sum(torque.values())
rate_req = np.degrees(A*w_bw) # deg/s

# ---- Servo linkage: servo degrees per gimbal degree (lever on the engine / servo horn) ----
tq_servo = d.get("servo_torque_Nm"); spd_servo = d.get("servo_speed_dps")
ratio_min = d.get("actuator_margin")*tq_req/tq_servo # below this the servo can't push hard enough
ratio_max = min(spd_servo/rate_req, d.get("servo_travel_max_deg")/rng) # above this it's too slow, or runs out of travel
ratio = d.get("gimbal_lever_r_m")/d.get("gimbal_horn_r_m") # set by the arm tab and the servo horn

print("-- Angle needed at hover (thrust = weight) --")
for case, (parts, tot, arm, T) in need.items():
    print(f"{case:8s} ({T/g:.1f} kg, arm {arm:.2f} m): " + ", ".join(f"{k} {v:.2f}" for k, v in parts.items()) + f" -> {tot:.2f} deg")
print(f"Gust side force:             {gust_F:.1f} N at {d.get('gust_mps'):.0f} m/s, {abs(x_cp - cases[worst][1]):.2f} m from the CG")
print(f"Worst case ({worst}):        {angle_req:.2f} deg needed, {use*rng:.1f} deg allowed ({use*100:.0f}% of +/-{rng:.0f} deg)")

print("\n-- What swings --")
print(f"Mass:                        {m_move:.2f} kg (engine {m_eng:.2f} + throttle valves, servos, fittings, igniter)")
print(f"Center of mass:              {(x_pivot - x_move)*1e3:.1f} mm below the pivot")
print(f"Inertia about the pivot:     {I_piv:.4f} kg m^2")

print("\n-- Gimbal axis (each of 2 axes) --")
for k, v in torque.items(): print(f"  {k:28s} {v:.3f} N m")
print(f"Torque needed:               {tq_req:.2f} N m; rate {rate_req:.0f} deg/s ({use*rng:.1f} deg sine at {d.get('actuator_bandwidth_hz'):.0f} Hz)")
print(f"Servo linkage window:        {ratio_min:.2f} to {ratio_max:.2f} servo deg per gimbal deg (torque x{d.get('actuator_margin'):.0f} / speed and travel)")
if True:
    print(f"Lever / horn ratio:          {ratio:.2f}: servo swings +/-{ratio*rng:.0f} deg for +/-{rng:.0f} deg; "
          f"gimbal gets {tq_servo*ratio:.1f} N m ({tq_servo*ratio/tq_req:.1f}x needed) and {spd_servo/ratio:.0f} deg/s ({spd_servo/ratio/rate_req:.1f}x needed)")

print("\n-- Swing clearance --")
r_body = d.get("body_id_m")/2 # body tube inside radius
r_nozzle = d.get("inlet_collar_outer_r_m") # widest part below the flange: the inlet collar
r_swing = r_nozzle + x_pivot*np.sin(rad(rng))
print(f"Nozzle end at +/-{rng:.0f} deg:       reaches r {r_swing*1e3:.0f} mm (collar r {r_nozzle*1e3:.0f} + {x_pivot*np.sin(rad(rng))*1e3:.0f} mm swing); "
      f"body inside r {r_body*1e3:.0f} mm")

# ---- Strength at full thrust x safety factor: two pins share the thrust at each joint ----
SF = d.get("gimbal_safety_factor"); Y = d.get("gimbal_yield_Pa"); P = F_max/2*SF
b_r = 0.5*(d.get("gimbal_ring_od_m") - d.get("gimbal_ring_id_m")); h_r = d.get("gimbal_ring_h_m"); d_b = d.get("gimbal_bush_od_m"); d_p = d.get("gimbal_pin_d_m")
M_ring = P*R_ring/2 # ring loaded at 4 alternating points: largest bending at the pins, where the bore weakens it
s_ring = M_ring*(h_r/2)/(b_r*(h_r**3 - d_b**3)/12)
t_ring = 0.207*P*R_ring/(0.25*b_r*h_r**2) # largest twist, halfway between the pins (0.25 ~ rectangular section)
s_pin = P/(np.pi/4*d_p**2) # shear across one pin
s_bush = P/(d_b*b_r) # bushing bearing on the ring
ry = d.get("gimbal_arm_riser_r_out_m") - d.get("gimbal_arm_riser_r_in_m"); rw = d.get("gimbal_arm_riser_w_m")
s_arm = P/(ry*(rw - d_p)) # riser in tension past the pin hole
t_wall = (d.get("body_od_m") - d.get("body_id_m"))/2; n_s = d.get("gimbal_screw_n")
s_tube = 2*P/(n_s*0.004*t_wall) # M4 screws bearing on the fiberglass tube wall
print(f"\n-- Strength at full thrust x {SF:.0f} (6061-T6 yield {Y/1e6:.0f} MPa) --")
print(f"Ring at the pins:            {s_ring/1e6:.0f} MPa bending, {t_ring/1e6:.0f} MPa twist")
print(f"Pin shear / bushing bearing: {s_pin/1e6:.0f} / {s_bush/1e6:.0f} MPa")
print(f"Arm riser tension:           {s_arm/1e6:.0f} MPa")
print(f"Body tube at the {n_s} screws:  {s_tube/1e6:.0f} MPa bearing")

# ---- Clearances at full tilt ----
th = rad(rng)
r_ring_o = d.get("gimbal_ring_od_m")/2; r_mount_i = d.get("gimbal_mount_id_m")/2
gap_mount = r_mount_i - (r_ring_o*np.cos(th) + h_r/2*np.sin(th))
r_arm_o = np.hypot(d.get("gimbal_arm_riser_r_out_m"), rw/2)
z_low_ = (h_r/2 + r_arm_o*np.sin(th))/np.cos(th) # the part of the riser below the pivot that swings level with the ring's lower inside edge
gap_arm = d.get("gimbal_ring_id_m")/2 - (r_arm_o*np.cos(th) + max(d.get("gimbal_arm_top_above_pivot_m"), z_low_)*np.sin(th))
# The flange sits below the ring, wider than its opening: check the gap as the engine swings on the inner pins
T_fl_ = d.get("inj_flange_thickness_m"); foot_ = d.get("gimbal_arm_foot_t_m"); r_bc_ = d.get("bolt_circle_d_m")/2
gap_flange = min((p - T_fl_)*np.cos(th) - d.get("flange_outer_d_m")/2*np.sin(th), # flange corner
                 (p - T_fl_ - foot_ - 0.005)*np.cos(th) - (r_bc_*np.cos(np.radians(22.5)) + 0.00425)*np.sin(th)) - h_r/2 # arm-foot bolt heads
window = (p - h_r/2 - R_ring*np.sin(th)) - d.get("inj_flange_thickness_m") # ethanol tube crossing above the flange, under the ring
tube_od = d.get("fu_line_od_m")
z_low = d.get("gimbal_bracket_t_m")/2 - d.get("servo_tab_above_bottom_m"); y_out = d.get("gimbal_servo_r_m") + d.get("servo_W_m")/2
r_servo = np.hypot(y_out, z_low) # outer corner of the servo case
print("\n-- Clearances at full tilt --")
print(f"Ring to mount:               {gap_mount*1e3:.1f} mm")
print(f"Arm to ring / flange to ring: {gap_arm*1e3:.1f} / {gap_flange*1e3:.1f} mm")
print(f"Ethanol tube under the ring: {window*1e3:.1f} mm of height for a {tube_od*1e3:.1f} mm tube")
print(f"Servo case corner:           r {r_servo*1e3:.1f} mm, body inside r {r_body*1e3:.1f} mm")
# Chamber outlet fitting: sticks out sideways from the port boss
x_fit = x_face - d.get("outlet_port_x_m") # height above the nozzle exit
r_fit = d.get("outlet_port_boss_top_r_m") + d.get("orb4_fitting_height_m")
r_fit_tilt = r_fit + (x_pivot - x_fit)*np.sin(th)
x_fit_top = x_fit + d.get("outlet_port_boss_d_m")/2 + r_fit*np.sin(th) # highest the fitting and its boss reach at full tilt
fit_in_tube = x_fit_top > d.get("body_tube_bottom_m") - 0.005 # within 5 mm of the tube's end
print(f"Outlet fitting at full tilt: r {r_fit_tilt*1e3:.0f} mm ({r_fit*1e3:.0f} mm straight)" + (f", body inside r {r_body*1e3:.1f} mm" if fit_in_tube else f", {(d.get('body_tube_bottom_m') - x_fit_top)*1e3:.0f} mm below the end of the body tube"))

problems = []
for name_, val in (("ring bending", s_ring), ("pin shear", s_pin), ("bushing bearing", s_bush), ("arm riser", s_arm)):
    if val > Y: problems.append(name_ + " over yield")
if t_ring > Y/np.sqrt(3): problems.append("ring twist over yield")
if s_tube > 100e6: problems.append("body tube bearing too high")
if min(gap_mount, gap_arm, gap_flange) < 0.0005: problems.append("gimbal parts touch at full tilt")
if window < tube_od: problems.append("ethanol tube doesn't fit under the ring")
if r_servo > r_body: problems.append("servo hits the body tube")
if fit_in_tube and r_fit_tilt > r_body - 0.002: problems.append("chamber outlet fitting hits the body tube")
if angle_req > use*rng: problems.append("gimbal range too small for the criteria")
if not ratio_min <= ratio <= ratio_max: problems.append("lever / horn ratio is outside the servo window")
if r_swing > r_body - 0.005: problems.append("engine swing hits the body")
print("\nProblems:                      " + ("; ".join(problems) if problems else "none"))

d.save_many({
    "gimbal_angle_needed_deg": (angle_req, "deg"),
    "gimbal_moving_mass_kg": (m_move, "kg"),
    "I_gimbal_pivot_kgm2": (I_piv, "kg*m^2"),
    "gimbal_torque_needed_Nm": (tq_req, "N*m"),
    "gimbal_rate_needed_dps": (rate_req, "deg/s"),
    "servo_ratio": (ratio, "-"),
}, "gimbal.py")
