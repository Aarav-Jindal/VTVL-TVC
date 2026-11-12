# Engine-off attitude check (rough), 2026-11-07: what the rocket's lean does between engine cutoff and the
# relight, with wind, a small lean and spin at cutoff, and an uncertain aerodynamic push. Then the gimbal's chance of taking it back.
# Also the other plan: the engine never goes off (minimum throttle after the climb), and the apogee that gives with the same propellant.
# 2-D pitch model: rigid body, gravity, axial drag, a slender-body push at whichever end leads (CNa, the unstable part), and crossflow
# drag along the body and folded legs (which also damps the swing). Masses and stations from structure.py, with the budget's margins.
# Rough on purpose: the coefficients are textbook estimates, so it runs every combination and reports the worst.
# Result (2026-11-07): it tumbles, so the flight plan was switched to the engine never off (flight_plan.py). Kept as the record.
import math
import numpy as np
import design as d
import structure as S
import flight_plan as fp

g, rho = 9.80665, 1.225
F, ve = d.get("thrust_N"), d.get("Isp_real_s")*g; mdot = F/ve
D = 2*S.R_O/1e3; A = math.pi*D**2/4
x_of = lambda u_mm: (u_mm + S.L_GAS)/1e3                                     # u (mm above the injector face) -> m above the nozzle exit
X_BASE, X_TOP = x_of(S.U_BOT), x_of(S.U_TOP); X_TIP = X_TOP + S.R_O/1e3
CD_UP, CD_DOWN = fp.CD_UP, fp.CD_DOWN                                        # axial drag: dome first (OpenRocket), base first (blunt end leading) [est]
ETA_CDC = 0.65*1.2                                                           # crossflow drag: cylinder 1.2 x finite-length factor 0.65 [Allen-Perkins, est]
CNA = (1.0, 2.0, 3.0)                                                        # slender-body push at the leading end, per rad: 2 by theory, spread for the blunt ends

# ---- The rocket at cutoff in the old plan (with margins): dry parts, helium, and the propellant kept for landing at the bottom of each tank ----
t_b = fp.engine_off_plan(d.get("m_wet_kg"), d.get("m_prop_kg"))[2]; burned = mdot*t_b
lumps = [(m, x_of(u)) for _, m, u in S.dry] + [(d.get("m_He_kg"), x_of(0.5*(S.PR.u_nb + S.PR.u_nt)))]
r_in = d.get("tank_outer_radius_m") - d.get("tank_wall_m"); A_in = math.pi*r_in**2
left = 1 - burned/d.get("m_prop_kg")
for tank, mk, rk in (("n2o_tank", "m_N2O_kg", "rho_N2O"), ("eth_tank", "m_ethanol_kg", "rho_ethanol")):
    h = d.get(mk)*left/(d.get(rk)*A_in); x0 = S.sec[tank][0] + S.cap/2       # sec is in m above the exit
    lumps.append((d.get(mk)*left, x0 + h/2))
m_cut = sum(m for m, _ in lumps); x_cg = sum(m*x for m, x in lumps)/m_cut
I_cut = sum(m*(x - x_cg)**2 for m, x in lumps)

# Body strips for crossflow: the tube (and dome) at full width, plus the folded legs (2 struts seen side-on) below
xs = np.arange(X_BASE, X_TIP, 0.05) + 0.025; widths = np.full_like(xs, D)
folded = S.LL.to_world(S.LL.CLOCKS[0], S.LL.A_FOLD, (S.LL.L_STRUT, 0, 0)); x_foot = x_of(-folded[0]); x_pin = x_of(S.LL.U_PIN)
xs_leg = np.arange(x_foot, x_pin, 0.05) + 0.025
xs = np.concatenate([xs, xs_leg]); widths = np.concatenate([widths, np.full_like(xs_leg, 2*0.030)])
s_strip = xs - x_cg; ds = 0.05
S_NOSE, S_BASE = (X_TOP + 0.05) - x_cg, X_BASE - x_cg                        # where the slender-body push acts: middle of the dome, or the blunt base

def climb(cd=CD_UP):
    """1-D climb at full thrust with drag: height and speed at cutoff."""
    h = v = 0.0; m = d.get("m_wet_kg"); dt = 0.002
    for _ in range(int(t_b/dt)):
        a = (F - 0.5*rho*cd*A*v*abs(v))/m - g; v += a*dt; h += v*dt; m -= mdot*dt
    return h, v

def deriv(st, wind, cna):
    X, Z, VX, VZ, th, w = st
    b = np.array([math.sin(th), math.cos(th)]); n = np.array([math.cos(th), -math.sin(th)])
    v_rel = np.array([wind - VX, -VZ])                                       # air relative to the CG
    a_ax = v_rel @ b; v_n = v_rel @ n - w*s_strip                            # crossflow at each strip, including the swing
    f_n = 0.5*rho*ETA_CDC*widths*ds*np.abs(v_n)*v_n                          # crossflow drag on each strip, along n
    s_lead = S_NOSE if a_ax < 0 else S_BASE                                  # air coming down the body: dome leads (going up); coming up it: base leads
    f_pot = 0.5*rho*A*cna*abs(a_ax)*(v_rel @ n)                              # slender-body push at the leading end, along n
    cd = CD_UP if a_ax < 0 else CD_DOWN
    f_ax = 0.5*rho*cd*A*abs(a_ax)*a_ax
    Fv = (f_n.sum() + f_pot)*n + f_ax*b
    M = (s_strip*f_n).sum() + s_lead*f_pot
    return np.array([VX, VZ, Fv[0]/m_cut, Fv[1]/m_cut - g, w, M/I_cut]), M

ARM = x_cg - d.get("x_gimbal_pivot_m"); TQ_GIMBAL = F*math.sin(math.radians(d.get("gimbal_range_deg")))*ARM
A_BRAKE = F/m_cut - g                                                        # full-thrust braking, upright, no drag (relight as late as allowed)

def fall(wind, th0, w0, cna, h0, v0, dt=0.002):
    st = np.array([0.0, h0, 0.0, v0, math.radians(th0), math.radians(w0)]); t = 0.0; th_max = 0.0; apogee = h0
    while True:
        k1, M = deriv(st, wind, cna); k2, _ = deriv(st + dt/2*k1, wind, cna); k3, _ = deriv(st + dt/2*k2, wind, cna); k4, _ = deriv(st + dt*k3, wind, cna)
        st = st + dt/6*(k1 + 2*k2 + 2*k3 + k4); t += dt
        apogee = max(apogee, st[1]); th_max = max(th_max, abs(math.degrees(st[4])))
        if st[3] < 0 and st[1] <= st[3]**2/(2*A_BRAKE): break                   # time to relight for a full-thrust landing burn
        if t > 60: break
    _, M = deriv(st, wind, cna)
    return dict(t_off=t, apogee=apogee, tilt=math.degrees(st[4]), rate=math.degrees(st[5]), th_max=th_max, h_relight=st[1],
                v_relight=-st[3], vx=st[2] - wind, M_aero=M)


if __name__ == "__main__":
    h0, v0 = climb()
    print("-- At engine cutoff (with margins) --")
    print(f"Mass {m_cut:.1f} kg, CG {x_cg:.3f} m above the nozzle exit, pitch inertia {I_cut:.0f} kg m^2. Cutoff at {t_b:.2f} s, {h0:.0f} m up, {v0:.1f} m/s")
    print(f"Gimbal at full thrust and full tilt: {TQ_GIMBAL:.0f} N m about the CG")
    print("\n-- Engine off, cutoff to relight (worst of CNa 1-3 and a lean of 0 or 2 deg and spin of 0 or 2 deg/s at cutoff, either way) --")
    print(f"{'wind m/s':>9s}{'engine off s':>13s}{'apogee m':>9s}{'relight m':>10s}{'fall m/s':>9s}{'worst lean':>11s}{'at relight':>11s}{'spin deg/s':>11s}{'aero vs gimbal':>16s}")
    rows = []
    for wind in (0.0, 4.0, 8.0, 12.0):
        runs = [fall(wind, th0, w0, cna, h0, v0) for cna in CNA for th0 in (0.0, 2.0, -2.0) for w0 in (0.0, 2.0, -2.0)]
        worst = max(runs, key=lambda r: abs(r["tilt"])); calm = fall(wind, 0.0, 0.0, 2.0, h0, v0)
        rows.append((wind, worst, calm))
        print(f"{wind:9.0f}{worst['t_off']:13.1f}{worst['apogee']:9.0f}{worst['h_relight']:10.0f}{worst['v_relight']:9.1f}{worst['th_max']:10.0f}d{abs(worst['tilt']):10.0f}d"
              f"{abs(worst['rate']):11.0f}{abs(worst['M_aero']):8.0f} / {TQ_GIMBAL:.0f}")
    print("(lean in deg from straight up; 'worst lean' is the most it reaches before the relight; 'aero' is the push turning it at the relight)")
    r = fp.always_lit_plan(d.get("m_wet_kg"), d.get("m_prop_kg"))
    print("\n-- The plan now (flight_plan.py): engine never off --")
    print(f"Full thrust {r['t1']:.1f} s, then {fp.THR_MIN*100:.0f}% over the top and down, then full thrust to land: apogee about {r['apogee']:.0f} m with drag "
          f"(same propellant, 5 s hover, 10% reserve). Flight {r['t_land']:.0f} s to the hover")
    d.save_many({"attitude_worst_lean_deg_8mps": (abs(rows[2][1]["tilt"]), "deg"), "attitude_worst_lean_deg_0mps": (abs(rows[0][1]["tilt"]), "deg"),
                 "always_lit_apogee_m": (r["apogee"], "m")}, "attitude_check.py")
