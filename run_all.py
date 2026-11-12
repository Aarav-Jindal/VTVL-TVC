# Runs every design script in the right order.
# cea_analysis.py runs twice: first for Isp (which the mass loop needs), then again with the new thrust
# from massIteration.py to size the engine. Isp doesn't depend on thrust, so one pass is enough.
import os
import subprocess
import sys

folder = os.path.dirname(os.path.abspath(__file__))
scripts = ["cea_analysis.py", "massIteration.py", "cea_analysis.py", "heat_transfer.py", "chamber_hardware.py", "injector.py", "injector_hardware.py", "feed_system.py", "vehicle_layout.py", "gimbal.py"]
for script in scripts:
    print(f"\n{'=' * 20} {script} {'=' * 20}", flush=True)
    result = subprocess.run([sys.executable, script], cwd=folder)
    if result.returncode != 0:
        print(f"\nStopped: {script} failed", flush=True)
        break
