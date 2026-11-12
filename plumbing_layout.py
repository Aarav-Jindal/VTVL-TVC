# Real stack heights from the real plumbing parts, for the CAD. Tanks keep the sizes the mass loop saved; the gaps come from what has to fit in them.
# Heights are u = mm above the injector face (engine frame x = -u). Run after massIteration.py. The mass budget doesn't see these lengths yet
# (the rocket is being built for 40-45 kg), so vehicle_layout.py's height is lower than this one.
import design as d
mm = lambda k: d.get(k)*1e3

GAP = 6.0          # visible tube between two tube-end parts
CAP = 35.0         # tube length each end cap uses (CAD)
BOSS_INSET = 8.0   # bottom cap's boss face, in from the tube end
FLEX_N2O = 165.0   # -10 hose between its ends: enough for its 165 mm bend radius at full tilt [Danfoss 2807-10]
HOSE_END = 55.0    # estimate: the hose end's length isn't published
BOTTLE_L = 515.6   # Luxfer ECLIPSE PRO L66P: 20.3 in long, 6.1 in across [maker]
GAP_BOTTLES = 180.0 # between the bottles: upper adapter, regulator, outlet adapter, connector, and the line turning out over the lower bottle's dome
parts = [("Connector, pintle -6 ORB", 39.1 - (mm("orb6_tap_depth_m") - 1)), ("gap", GAP), ("OTV SS-1RS8", 71.1), ("gap", GAP),
         ("Hose end", HOSE_END), ("-10 hose flex", FLEX_N2O), ("Hose end", HOSE_END), ("gap", GAP), ("MOV SS-43GS8", 71.1), ("gap", GAP), ("F-1 SS-8F-40", 88.6), ("gap", GAP), ("Connector, cap -8 ORB", 41.9 - 13.3)]
u = mm("inj_hub_thickness_m") + mm("pintle_head_len_m") # top of the pintle
stations = {}
for name, length in parts:
    if name != "gap": stations[name] = (u, u + length)
    u += length
u_n2o_bot = u - BOSS_INSET
L_n2o = mm("L_tank_N2O_m") + 2*CAP; L_eth = mm("L_tank_ethanol_m") + 2*CAP
GAP_N2O_ETH = 230.0 # ethanol filter and main valve on the axis, nitrous vent and burst items around them
GAP_ETH_HE = 270.0  # adapter, regulator (size estimated), tee, over the ethanol top cap fittings
u_n2o_top = u_n2o_bot + L_n2o
u_eth_bot = u_n2o_top + GAP_N2O_ETH; u_eth_top = u_eth_bot + L_eth
u_he_port = u_eth_top + GAP_ETH_HE; u_he2_port = u_he_port + BOTTLE_L + GAP_BOTTLES; u_he_end = u_he2_port + BOTTLE_L
u_top = u_he_end + mm("section_gap_m") + mm("avionics_bay_len_m")
L_gas = d.get("L_gas_total_mm")
x_tube_end = L_gas - d.get("body_tube_bottom_m")*1e3 # engine frame
print("-- Real stack (mm above the injector face) --")
for k, v in stations.items(): print(f"{k:28s} {v[0]:7.1f} to {v[1]:7.1f}")
print(f"{'Nitrous tank tube':28s} {u_n2o_bot:7.1f} to {u_n2o_top:7.1f}  ({L_n2o:.0f} mm)")
print(f"{'Ethanol tank tube':28s} {u_eth_bot:7.1f} to {u_eth_top:7.1f}  ({L_eth:.0f} mm)")
print(f"{'Helium bottles':28s} {u_he_port:7.1f} to {u_he_port + BOTTLE_L:7.1f}, {u_he2_port:7.1f} to {u_he_end:7.1f}")
print(f"{'Top of the body':28s} {u_top:7.1f}")
L_body = u_top + x_tube_end
print(f"Height above the nozzle exit: {(u_top + L_gas)/1e3:.2f} m (mass budget's layout: {d.get('L_vehicle_m'):.2f} m). Body tube {L_body/1e3:.2f} m = {L_body/25.4:.0f} in, in lengths of 60 in at most: {int(L_body/25.4//60) + 1} sections")
d.save_many({"plumb_u_n2o_bot_mm": (u_n2o_bot, "mm"), "plumb_u_n2o_top_mm": (u_n2o_top, "mm"), "plumb_u_eth_bot_mm": (u_eth_bot, "mm"), "plumb_u_eth_top_mm": (u_eth_top, "mm"),
             "plumb_u_he_port_mm": (u_he_port, "mm"), "plumb_u_he2_port_mm": (u_he2_port, "mm"), "plumb_bottle_len_mm": (BOTTLE_L, "mm"), "plumb_u_top_mm": (u_top, "mm"), "plumb_L_n2o_tube_mm": (L_n2o, "mm"), "plumb_L_eth_tube_mm": (L_eth, "mm"),
             "plumb_u_mov_mid_mm": (sum(stations["MOV SS-43GS8"])/2, "mm"), "plumb_u_f1_mid_mm": (sum(stations["F-1 SS-8F-40"])/2, "mm"), "plumb_u_otv_mid_mm": (sum(stations["OTV SS-1RS8"])/2, "mm"),
             }, "plumbing_layout.py")
