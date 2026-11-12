"""
Shared design values for the whole project.

Every script reads its inputs from design_values.json and writes its results back to it,
so there is one place that always holds the current design.

    import design as d
    Pc = d.get("Pc_psi")                                          # read one value
    d.save("R_throat_mm", R_throat*1e3, "mm", "cea_analysis.py")  # write one value
    d.save_many({"Isp_real_s": (Isp_real, "s"),                   # write several at once
                 "cstar_real_mps": (cstar_real, "m/s")}, "cea_analysis.py")

Inputs have source "input": those are the values you choose, edit them in design_values.json.
Results have source = the script that calculates them, and get overwritten every time it runs.
"""
import json
import os

FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "design_values.json")


def load():
    with open(FILE, encoding="utf-8") as f:
        return json.load(f)


def get(name):
    data = load()
    if name not in data:
        raise KeyError(f"'{name}' is not in design_values.json yet. Run the script that calculates it first.")
    return data[name]["value"]


def save_many(values, source):
    data = load()
    for name, (value, unit) in values.items():
        if hasattr(value, "item"): # numpy number -> plain Python number
            value = value.item()
        entry = {"value": value, "unit": unit, "source": source}
        if "note" in data.get(name, {}): # keep any note written next to the value
            entry["note"] = data[name]["note"]
        data[name] = entry
    with open(FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def save(name, value, unit, source):
    save_many({name: (value, unit)}, source)
