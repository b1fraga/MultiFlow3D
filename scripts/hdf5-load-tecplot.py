import glob
import os
import re

import tecplot as tp
from tecplot.constant import *

tp.session.connect()

data_dir = os.path.dirname(os.path.abspath(__file__))

# Find all files for zone 0001
files = sorted(glob.glob(os.path.join(data_dir, "tecout_*.h5")))

dataset = None

for f in files:
    print(f"Processing: {f}")
    
    # Extract timestep from filename, e.g. tecout_0001_0003.h5 -> 3
    m = re.search(r'tecout_(\d+)_\d+\.h5$', os.path.basename(f))
    timestep = int(m.group(1))

    tp.macro.execute_command(f"""$!ReadDataSet '"-F" "1" "{f}" "-D" "4" "W" "P" "U" "V" "-R" "x" "y" "z" "-K" "1" "1" "1"'
      DataSetReader = 'HDF5 Loader'
      ReadDataOption = {'New' if dataset is None else 'Append'}
      ResetStyle = No
      AssignStrandIDs = Yes
      InitialPlotType = Cartesian3D
      InitialPlotFirstZoneOnly = No
      AddZonesToExistingStrands = No
      VarLoadMode = ByName""")

    dataset = tp.active_frame().dataset

    zone = dataset.zone(dataset.num_zones - 1)
    zone.name = f"zone {dataset.num_zones-1}"

    zone.strand = 1
    zone.solution_time = timestep

tp.active_frame().plot().fieldmaps(0).show=True
tp.macro.execute_command('$!GlobalTime SolutionTime = 1')
tp.active_frame().plot().view.fit(consider_blanking=True)
