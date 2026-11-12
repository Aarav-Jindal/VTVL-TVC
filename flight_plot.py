# Whole-flight graph for the flight plan (flight_plan.py: engine never off): height, speed and throttle against time, with and without
# the mass margins. OpenRocket only covers the climb (it can't steer a powered descent), so this 1-D plan is the picture of the rest.
# Saves sim/flight_profile.png.
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import design as d
import structure as S
import flight_plan as fp

HERE = os.path.dirname(os.path.abspath(__file__))
m_prop, m_He = d.get("m_prop_kg"), d.get("m_He_kg")
cases = [("With margins", d.get("m_wet_kg"), "#2a6f97"),
         ("No margins", d.get("m_wet_kg") - d.get("m_dry_kg")*(1 - 1/S.K_MARGIN), "#c8553d")]
fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
for label, wet, color in cases:
    plan = fp.always_lit(wet, fp.always_lit_plan(wet, m_prop)["t1"], record=True)
    t, h, v, thr, m = (np.array(c) for c in zip(*plan["rows"]))
    k = slice(None, None, 50)
    ax[0].plot(t[k], h[k], color=color, label=f"{label} ({wet:.1f} kg): apogee {plan['apogee']:.0f} m")
    ax[1].plot(t[k], v[k], color=color)
    ax[2].plot(t[k], thr[k]*100, color=color)
    for a in ax: a.axvline(plan["t_land"], color=color, lw=0.8, ls=":")
ax[0].set_ylabel("Height (m)"); ax[1].set_ylabel("Speed up (m/s)"); ax[2].set_ylabel("Throttle (%)"); ax[2].set_xlabel("Time (s)")
ax[1].axhline(0, color="0.6", lw=0.8); ax[2].set_ylim(0, 110)
ax[0].legend(loc="upper right", frameon=False)
ax[0].set_title("Flight plan: engine never off (full thrust, 40% over the top and down, full thrust to land, 5 s hover)\n"
                "1-D, Cd 0.52 going up and 0.8 coming down; dotted lines: landing burn done (stopped about 10 m up), then the hover", fontsize=10)
for a in ax: a.grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(HERE, "sim", "flight_profile.png"); fig.savefig(out, dpi=130)
print("saved", out)
