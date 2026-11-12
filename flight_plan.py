# The flight plan: how the propellant is spent, shared by openrocket_sim.py, flight_plot.py and attitude_check.py.
# Current plan: the engine never goes off. Full thrust to climb, then the minimum throttle over the top
# and on the way down, then full thrust to stop at the ground, then a 5 s hover, with 10% of the propellant kept in reserve.
# The gimbal steers the whole flight; attitude_check.py shows that a plan with the engine off from cutoff to a relight tumbles.
# That plan is engine_off_plan(), used by attitude_check.py. 1-D, straight up and down.
import math
import design as d

g, RHO = 9.80665, 1.225
F = d.get("thrust_N"); VE = d.get("Isp_real_s")*g
THR_MIN = d.get("throttle_min"); HOVER_S = 5.0; RESERVE = 1.1
A_REF = math.pi*(d.get("body_od_m")/2)**2
CD_UP, CD_DOWN = 0.52, 0.8          # dome first (OpenRocket); base first, blunt end leading [est]

def always_lit(m0, t1, drag=True, dt=0.002, record=False, prop=None):
    """Fly the plan with t1 s of full thrust. Returns apogee, propellant used (with the hover), phase times and, if asked, the timeline.
    prop: stop counting once this much is gone (a too-long climb can leave the rocket light enough to rise on minimum throttle forever)."""
    m, h, v, t, apo, t_apo, phase = m0, 0.0, 0.0, 0.0, 0.0, 0.0, 1
    t_brake = None; rows = []
    while True:
        if prop is not None and m0 - m > prop:                               # ran dry: report it as too much
            return dict(apogee=apo, t_apogee=t_apo, used=float("inf"), t1=t1, t_brake=None, t_land=None, t_end=None, rows=rows)
        T = F if phase in (1, 3) else THR_MIN*F
        D = 0.5*RHO*(CD_UP if v >= 0 else CD_DOWN)*A_REF*v*abs(v) if drag else 0.0
        if record: rows.append((t, h, v, T/F, m))
        v += ((T - D)/m - g)*dt; h += v*dt; m -= T/VE*dt; t += dt
        if h > apo: apo, t_apo = h, t
        if phase == 1 and t >= t1: phase = 2
        if phase == 2 and v < 0 and h <= v*v/(2*(F/m - g)): phase, t_brake = 3, t   # brake at full thrust to a stop at the ground
        if phase == 3 and v >= 0: break
    t_land = t
    if record:                                                               # the hover: thrust = weight, falling as the mass goes
        for _ in range(int(HOVER_S/dt)):
            rows.append((t, h, 0.0, m*g/F, m)); m -= m*g/VE*dt; t += dt
    m_end = m if record else m*math.exp(-g*HOVER_S/VE)
    return dict(apogee=apo, t_apogee=t_apo, used=m0 - m_end, t1=t1, t_brake=t_brake, t_land=t_land, t_end=t_land + HOVER_S, rows=rows)

def always_lit_plan(m0, prop, drag=True):
    """The longest full-thrust climb the propellant allows, keeping the reserve."""
    lo, hi = 0.5, 40.0; best = None
    for _ in range(40):
        t1 = 0.5*(lo + hi); r = always_lit(m0, t1, drag, prop=prop)
        if r["used"] > prop/RESERVE: hi = t1
        else: lo, best = t1, r
    return best

def engine_off_plan(m0, prop):
    """Engine-off plan, no drag: full thrust to a burnout speed, engine off over the top and
    down, relit for a full-thrust landing burn, 5 s hover, 10% reserve. Returns burnout speed, apogee, burn time, propellant burned."""
    mdot = F/VE
    def flight(vb):
        dt = 0.001; m = m0; v = h = t = 0.0
        while v < vb: v += (F/m - g)*dt; h += v*dt; m -= mdot*dt; t += dt
        ap = h + v*v/(2*g); lo, hi = 0.0, 250.0
        for _ in range(30):
            vi = 0.5*(lo + hi); x = ap - vi*vi/(2*g); mm_ = m; vv = vi
            while vv > 0: vv -= (F/mm_ - g)*0.01; x -= vv*0.01; mm_ -= mdot*0.01
            lo, hi = (vi, hi) if x > 0 else (lo, vi)
        return ap, m0 - mm_*math.exp(-g*HOVER_S/VE), t, m0 - m
    lo, hi = 1.0, 250.0; best = None
    for _ in range(30):
        vb = 0.5*(lo + hi); ap, used, t, burned = flight(vb)
        if used > prop/RESERVE: hi = vb
        else: lo, best = vb, (vb, ap, t, burned)
    return best

if __name__ == "__main__":
    m0, prop = d.get("m_wet_kg"), d.get("m_prop_kg")
    for drag in (False, True):
        r = always_lit_plan(m0, prop, drag)
        print(f"{'with' if drag else 'no  '} drag: full thrust {r['t1']:.2f} s, {THR_MIN*100:.0f}% from there, apogee {r['apogee']:.0f} m at {r['t_apogee']:.1f} s, "
              f"braking at {r['t_brake']:.1f} s, landed {r['t_land']:.1f} s, hover ends {r['t_end']:.1f} s, propellant {r['used']:.2f} of {prop:.2f} kg")
