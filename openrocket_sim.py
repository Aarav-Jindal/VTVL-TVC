# OpenRocket model of the climb: builds the rocket from the same mass layout as structure.py (vehicle_layout.py's stations,
# the CAD's tube and dome), writes the flight's thrust as a motor, runs the simulations and saves .ork files in sim/.
# Runs the dome nose, plus ogive and von Karman noses for comparison (sim/nose_variants).
# The flight plan is flight_plan.py's (engine never off): full thrust, then minimum throttle over the top and
# down, then full thrust to land, then the hover. The motor carries that whole thrust history, so the mass drops as it should, but the
# runs stop a few seconds past apogee: OpenRocket can't steer, and a rocket that's stable going up is unstable coming down tail first.
# Zero wind, launched straight up; the gimbal is what holds it straight in reality, so wind runs would mean nothing here.
# Needs OpenRocket 24.12 (its own Java) and JPype + numpy in the Python that runs this (pip install jpype1 numpy).
# Writes its motors into OpenRocket's user motor folder (%APPDATA%/OpenRocket/ThrustCurves) so the saved files open with them.
import glob, math, os, shutil
import numpy as np
import design as d
import structure as S                                                       # the mass layout, stations and margins (prints nothing on import)
import flight_plan as fp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sim"); os.makedirs(os.path.join(OUT, "nose_variants"), exist_ok=True)
MOTOR_DIR = os.path.join(os.environ["APPDATA"], "OpenRocket", "ThrustCurves")
g0 = 9.80665
F, Isp = d.get("thrust_N"), d.get("Isp_real_s"); ve = Isp*g0; mdot = F/ve
m_prop, m_He = d.get("m_prop_kg"), d.get("m_He_kg")
U_TOP, U_BOT = S.U_TOP, S.U_BOT
R_O, R_I = S.R_O, S.R_I                                                    # mm

# ---- Mass groups (kg, u mm), no margin. structure.py's lumps regrouped by where the parts sit, so OpenRocket's parts list reads sensibly ----
rows, kg, spots, u_of = S.rows, S.kg, S.spots, S.u_of
NAMES = {"engine_top": "Throttle valves, actuators, torch, PT-8", "engine": "Engine bolts and port fittings", "gimbal": "Gimbal hardware and servos",
         "engine_bay": "Engine bay: hoses, spark module, low BEC", "legs": "Landing legs, folded", "n2o_bottom": "Nitrous tank bottom: MOV, filter, fill",
         "n2o_top": "Nitrous tank top: vent, relief, sensors", "eth_bottom": "Ethanol tank bottom: MFV, filter, fill", "eth_top": "Ethanol tank top: vent, relief, sensors",
         "he_valves": "Helium regulators and valves", "he_bottle": "Bottle rings", "avionics": "Avionics bay"}
groups = {}
def add(name, m, u):
    m0, mu = groups.get(name, (0.0, 0.0)); groups[name] = (m0 + m, mu + m*u)
spread = d.get("m_body_tube_kg")
for r in rows:
    w, it = r["where"], r["item"]
    if it.startswith("Dome"): add("Dome inserts and screws", 0.020, S.AV.U_DOME_SCREW); continue  # the shell is OpenRocket's nose cone (PETG, 2 mm)
    if w == "body" and it.startswith("Body tube joints"):
        for u in S.AV.U_JOINT: add(f"Body joint u {u:.0f}", kg(r)/2, u)
    elif w == "body": spread += kg(r)
    elif w == "rings":
        for u in S.ML.U_TANK_RINGS: add(f"Tank ring u {u:.0f}", kg(r)/4, u)
    else: add(NAMES[w], kg(r), u_of(spots[w]))
add("Engine (CAD)", d.get("m_engine_cad_kg"), u_of(spots["engine"])); add("Gimbal (CAD)", d.get("m_gimbal_cad_kg"), S.U_PIVOT)
add("Nitrous tank (tube + caps)", d.get("m_tank_N2O_kg"), 0.5*(S.PR.u_nb + S.PR.u_nt)); add("Ethanol tank (tube + caps)", d.get("m_tank_ethanol_kg"), 0.5*(S.PR.u_eb + S.PR.u_et))
for i, u in enumerate(S.U_BOTTLES): add(f"Helium bottle {('lower', 'upper')[i]}", S.m_bottle, u)
for u in S.U_BOTTLES: add("Helium gas", m_He/2, u)                          # moves to the tanks during the flight; ignored here
DRY_NO_MARGIN = sum(m for n, (m, _) in groups.items() if n != "Helium gas") + spread + 0.212  # 0.212: OpenRocket's PETG dome, about the CAD's 0.218

# Propellant columns: liquid from each tank's bottom up (vehicle_layout.py). The climb burns the top of each column; the rest stays as fixed mass
r_in = d.get("tank_outer_radius_m") - d.get("tank_wall_m"); A_in = math.pi*r_in**2
def column(tank, m, rho):
    return u_of(S.sec[tank][0] + S.cap/2), m/(rho*A_in)*1e3              # bottom u, height mm
COL = {"Nitrous": (column("n2o_tank", d.get("m_N2O_kg"), d.get("rho_N2O")), d.get("m_N2O_kg")),
       "Ethanol": (column("eth_tank", d.get("m_ethanol_kg"), d.get("rho_ethanol")), d.get("m_ethanol_kg"))}

# ---- OpenRocket ----
def write_motor(name, plan):
    """The whole flight's thrust as one motor: full, minimum throttle, full to land, then the hover (thrust = weight, every 0.5 s)."""
    t1, tb, tl = plan["t1"], plan["t_brake"], plan["t_land"]; e = 0.02
    pts = [(e, F), (t1, F), (t1 + e, fp.THR_MIN*F), (tb, fp.THR_MIN*F), (tb + e, F), (tl, F)]
    hover = [(t, T*F) for t, h, v, T, m in plan["rows"] if t > tl + e]
    pts += hover[::250] + [(plan["t_end"], hover[-1][1]), (plan["t_end"] + e, 0.0)]
    path = os.path.join(MOTOR_DIR, f"{name}.eng")
    with open(path, "w") as f:
        f.write(f"; VTVL flight (flight_plan.py): {F:.0f} N for {t1:.2f} s, {fp.THR_MIN*100:.0f}% to {tb:.2f} s, full to land at {tl:.2f} s, "
                f"{fp.HOVER_S:.0f} s hover; {plan['used']:.3f} kg burned (openrocket_sim.py)\n")
        f.write(f"{name} 150 500 P {plan['used']:.4f} {plan['used'] + 0.001:.4f} VTVL\n")
        for t, T in pts: f.write(f"   {t:.3f} {T:.1f}\n")
        f.write(";\n")
    shutil.copy(path, OUT)
    return path

for old in glob.glob(os.path.join(MOTOR_DIR, "VTVL-climb-*.eng")) + glob.glob(os.path.join(OUT, "VTVL-climb-*.eng")):
    os.remove(old)                                                           # leftover motor files from the engine-off plan
CASES = {}
for key, label, K in (("margin", "with margins", S.K_MARGIN), ("no_margin", "no margins", 1.0)):
    wet = DRY_NO_MARGIN*K + m_prop + m_He
    plan = fp.always_lit(wet, fp.always_lit_plan(wet, m_prop)["t1"], record=True)
    nd = fp.always_lit_plan(wet, m_prop, drag=False)
    CASES[key] = dict(label=label, K=K, wet=wet, plan=plan, ap0=nd["apogee"], burned=plan["used"], motor=f"VTVL-flight-{key.replace('_', '')}")
    write_motor(CASES[key]["motor"], plan)

import jpype, jpype.imports                                                 # after the motor files exist, so OpenRocket's loader finds them
jpype.startJVM(r"C:\Program Files\OpenRocket\jre\bin\server\jvm.dll", classpath=[r"C:\Program Files\OpenRocket\jar\OpenRocket-24.12.jar"], convertStrings=False)
from info.openrocket.swing.startup import GuiModule
from info.openrocket.core.plugin import PluginModule
from info.openrocket.core.startup import Application
from com.google.inject import Guice
_gui = GuiModule(); Application.setInjector(Guice.createInjector(_gui, PluginModule())); _gui.startLoader()
from info.openrocket.core.rocketcomponent import Rocket, AxialStage, NoseCone, Transition, BodyTube, InnerTube, MassComponent, PodSet, FlightConfigurationId, ExternalComponent, TrapezoidFinSet
from info.openrocket.core.rocketcomponent.position import AxialMethod, RadiusMethod
from info.openrocket.core.document import OpenRocketDocumentFactory, Simulation
from info.openrocket.core.simulation import FlightDataType
from info.openrocket.core.motor import MotorConfiguration
from info.openrocket.core.material import Material
from info.openrocket.core.file import GeneralRocketSaver
from info.openrocket.core.simulation.extension.example import StopSimulation
from java.io import File
T_STOP = int(math.ceil(max(c["plan"]["t_apogee"] for c in CASES.values()) + 3)) # s: a few seconds past apogee (see the top)

def find_motor(name):
    for ms in Application.getMotorSetDatabase().getMotorSets():
        for mo in ms.getMotors():
            if str(mo.getDesignation()) == name: return mo
    raise RuntimeError(f"motor {name} not found in OpenRocket's database")
FINISH = next(f for f in ExternalComponent.Finish.values() if "Regular paint" in str(f))
def bulk(name, rho): return Material.newMaterial(Material.Type.BULK, name, float(rho), True)

def build(case, shape, ratio, param=0.0, stand_in=True):
    """The rocket for one mass case and nose. ratio = nose length / body diameter. stand_in: the gimbal stand-in fins (see below)."""
    c = CASES[case]; K = c["K"]
    rocket = Rocket(); rocket.setName(f"VTVL, {c['label']}")
    stage = AxialStage(); stage.setName("Rocket"); rocket.addChild(stage)
    nose = NoseCone(); nose.setName("Nose"); nose.setShapeType(shape); nose.setShapeParameter(float(param))
    nose.setLength(ratio*2*R_O/1e3); nose.setBaseRadius(R_O/1e3); nose.setThickness(0.002)
    nose.setAftShoulderLength(0.030); nose.setAftShoulderRadius(R_I/1e3); nose.setAftShoulderThickness(0.002)
    nose.setMaterial(bulk("PETG (printed)", 1270*K)); nose.setFinish(FINISH); stage.addChild(nose)
    tube = BodyTube(); tube.setName("Body tube (G12) + spread parts"); tube.setLength((U_TOP - U_BOT)/1e3)
    tube.setOuterRadius(R_O/1e3); tube.setThickness((R_O - R_I)/1e3); tube.setFinish(FINISH)
    tube.setMassOverridden(True); tube.setOverrideMass(spread*K); stage.addChild(tube)
    def mass(name, m, u, L=0.05):
        mc = MassComponent(); mc.setName(name); mc.setLength(L); mc.setRadius(0.06); mc.setComponentMass(float(m)); tube.addChild(mc)
        mc.setAxialMethod(AxialMethod.TOP); mc.setAxialOffset((U_TOP - u)/1e3 - L/2)  # OpenRocket places a part only once it has a parent
    for n, (m, mu) in groups.items():
        mass(n, m*(1 if n == "Helium gas" else K), mu/m)
    left = m_prop - c["burned"]; u_burn = 0.0
    for liq, ((u0, h), m_l) in COL.items():                                # the 10% reserve stays at the bottom of each column
        f_left = left/m_prop; h_left = h*f_left
        mass(f"{liq} reserve", m_l*f_left, u0 + h_left/2)
        u_burn += (m_l*(1 - f_left))*(u0 + (h_left + h)/2)
    u_burn /= c["burned"]
    mount = InnerTube(); mount.setName("Propellant burned in flight (motor)"); mount.setLength(0.5); mount.setOuterRadius(0.0762); mount.setThickness(0.0004)
    tube.addChild(mount); mount.setAxialMethod(AxialMethod.TOP); mount.setAxialOffset((U_TOP - u_burn)/1e3 - 0.25)
    mount.setMassOverridden(True); mount.setOverrideMass(0.0); mount.setMotorMount(True)
    pods = PodSet(); pods.setName("Folded legs (drag only)"); tube.addChild(pods); pods.setInstanceCount(4); pods.setRadiusMethod(RadiusMethod.FREE)
    pods.setRadiusOffset(0.182); pods.setAxialMethod(AxialMethod.TOP); pods.setAxialOffset((U_TOP - S.LL.U_PIN)/1e3)
    cap = NoseCone(); cap.setName("Hinge end"); pods.addChild(cap); cap.setShapeType(Transition.Shape.ELLIPSOID); cap.setLength(0.015)
    cap.setBaseRadius(0.015); cap.setThickness(0.0015); cap.setMassOverridden(True); cap.setOverrideMass(0.0)
    strut = BodyTube(); strut.setName("Strut"); pods.addChild(strut); strut.setLength(S.LL.L_STRUT/1e3); strut.setOuterRadius(0.015); strut.setThickness(0.0015)
    strut.setFinish(FINISH); strut.setMassOverridden(True); strut.setOverrideMass(0.0)
    if stand_in:
        # With no fins the rocket is unstable (center of pressure far ahead of the CG), and OpenRocket can't model the gimbal that holds
        # it straight: unchecked, it tumbles about 10 s in. These fins only stand in for the gimbal: drag and mass set to zero, and flying
        # straight up with no wind they never see an angle, so they make no force. They are not a design.
        fins = TrapezoidFinSet(); tube.addChild(fins); fins.setName("Gimbal stand-in (no drag, no mass, not real)")
        fins.setFinCount(4); fins.setRootChord(0.30); fins.setTipChord(0.15); fins.setSweep(0.15); fins.setHeight(0.30); fins.setThickness(0.003)
        fins.setAxialMethod(AxialMethod.BOTTOM); fins.setAxialOffset(0.0); fins.setAngleOffset(math.radians(45))
        fins.setMassOverridden(True); fins.setOverrideMass(0.0); fins.setCDOverridden(True); fins.setOverrideCD(0.0)
    fcid = FlightConfigurationId(); rocket.createFlightConfiguration(fcid)
    mc = MotorConfiguration(mount, fcid); mc.setMotor(find_motor(c["motor"])); mount.setMotorConfig(mc, fcid)
    rocket.setSelectedConfiguration(fcid)
    return rocket, fcid

def simulate(rocket, fcid, name, t_max=None):
    doc = OpenRocketDocumentFactory.createDocumentFromRocket(rocket)
    sim = Simulation(doc, rocket); sim.setName(name); sim.setFlightConfigurationId(fcid)
    o = sim.getOptions(); o.setLaunchRodLength(1.0); o.setLaunchRodAngle(0.0); o.setWindSpeedAverage(0.0); o.setWindSpeedDeviation(0.0); o.setWindTurbulenceIntensity(0.0)
    stop = StopSimulation(); stop.setStopTime(int(t_max or T_STOP)); sim.getSimulationExtensions().add(stop) # saved with the file, so reruns stop too
    sim.simulate(); doc.addSimulation(sim)
    fd = sim.getSimulatedData(); b = fd.getBranch(0)
    return doc, fd, (lambda t: np.array([float(v) for v in b.get(t)]))

def run(case, shape, ratio, param, label, save_as):
    _, _, colb = simulate(*build(case, shape, ratio, param, stand_in=False), "bare", t_max=3.0) # the rocket's own stability, 2 s in
    tb_ = colb(FlightDataType.TYPE_TIME); bare = colb(FlightDataType.TYPE_STABILITY)[int(np.argmin(abs(tb_ - 2.0)))]
    rocket, fcid = build(case, shape, ratio, param)
    doc, fd, col = simulate(rocket, fcid, f"Climb, {CASES[case]['label']}, no wind")
    t, v, mass_, cd, aoa, drag = (col(getattr(FlightDataType, k)) for k in ("TYPE_TIME", "TYPE_VELOCITY_Z", "TYPE_MASS", "TYPE_DRAG_COEFF", "TYPE_AOA", "TYPE_DRAG_FORCE"))
    tb = CASES[case]["plan"]["t1"]; ib = int(np.argmin(abs(t - (tb - 0.1)))); ic = int(np.argmin(abs(t - (tb + 1.0))))
    t_ap = float(fd.getTimeToApogee()); up = t <= t_ap; burn = (t > 1.0) & (t <= t_ap - 2.0)
    res = dict(label=label, case=case, apogee=float(fd.getMaxAltitude()), vmax=float(fd.getMaxVelocity()), amax=float(fd.getMaxAcceleration()),
               t_ap=t_ap, m0=mass_[0], v_b=v[ib], cd_burn=cd[ib], cd_coast=cd[ic], drag_max=np.nanmax(drag[up]), aoa_max=np.degrees(np.nanmax(np.abs(aoa[burn]))),
               stab=bare, nose_kg=float(rocket.getChild(0).getChild(0).getComponentMass()), warnings=str(fd.getWarningSet()))
    GeneralRocketSaver().save(File(save_as), doc)
    return res

if __name__ == "__main__":
    print("-- Flight plan (flight_plan.py: engine never off; 10% of the propellant kept in reserve) --")
    for c in CASES.values():
        p = c["plan"]
        print(f"{c['label']:13s} wet {c['wet']:.2f} kg: full thrust {p['t1']:.2f} s, {fp.THR_MIN*100:.0f}% to {p['t_brake']:.1f} s, full to land at {p['t_land']:.1f} s, "
              f"hover to {p['t_end']:.1f} s; {c['burned']:.2f} of {m_prop:.2f} kg; apogee {p['apogee']:.0f} m (1-D, Cd {fp.CD_UP}), {c['ap0']:.0f} m with no drag")
    Shape = Transition.Shape
    noses = [("Dome (now)", Shape.ELLIPSOID, 0.5, 0.0), ("Ogive 2:1", Shape.OGIVE, 2.0, 1.0), ("Ogive 3:1", Shape.OGIVE, 3.0, 1.0),
             ("Von Karman 2:1", Shape.HAACK, 2.0, 0.0), ("Von Karman 3:1", Shape.HAACK, 3.0, 0.0)]
    results = []
    for case in CASES:
        for label, shape, ratio, param in noses:
            fname = "VTVL_dome.ork" if (case == "margin" and label.startswith("Dome")) else \
                    os.path.join("nose_variants", f"VTVL_{label.split(' (')[0].replace(' ', '_').replace(':', 'to')}_{case}.ork")
            results.append(run(case, shape, ratio, param, label, os.path.join(OUT, fname)))
    print("\n-- OpenRocket, no wind, straight up --")
    print(f"{'nose':16s}{'mass case':14s}{'nose kg':>8s}{'liftoff kg':>11s}{'apogee m':>10s}{'max m/s':>9s}{'apogee s':>10s}{'Cd full':>9s}{'Cd 40%':>10s}{'max drag N':>11s}{'max AoA':>9s}{'bare stab cal':>15s}")
    for r in results:
        print(f"{r['label']:16s}{CASES[r['case']]['label']:14s}{r['nose_kg']:8.3f}{r['m0']:11.2f}{r['apogee']:10.0f}{r['vmax']:9.1f}{r['t_ap']:10.1f}{r['cd_burn']:9.3f}{r['cd_coast']:10.3f}{r['drag_max']:11.1f}{r['aoa_max']:9.3f}{r['stab']:15.2f}")
    print("(Drag losses: no-drag apogee " + ", ".join(f"{c['label']} {c['ap0']:.0f} m" for c in CASES.values()) + ". Bare stability: without the stand-in fins, 2 s in; negative = unstable)")
    for r in results:
        if r["warnings"] not in ("[]", ""): print(f"Warnings, {r['label']}, {CASES[r['case']]['label']}: {r['warnings']}")
    dome = next(r for r in results if r["case"] == "margin" and r["label"].startswith("Dome"))
    d.save_many({"or_apogee_margin_m": (dome["apogee"], "m"), "or_apogee_no_margin_m": (next(r["apogee"] for r in results if r["case"] == "no_margin" and r["label"].startswith("Dome")), "m"),
                 "or_stability_cal": (dome["stab"], "cal"), "or_full_thrust_s_margin": (CASES["margin"]["plan"]["t1"], "s")}, "openrocket_sim.py")
    jpype.shutdownJVM()
