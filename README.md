# py-1-DT
plz add some te\xt
write sometihings
change
## Digital Twin Aquarium Dashboard

This Python application builds a reduced-order model from ANSYS simulation exports and displays predicted temperature, velocity, and dissolved oxygen in a 3D PyQt dashboard. It also supports probe readings, alarm indicators, experimental validation, and exporting validation results to Excel.

### Setup

Install Python and the libraries used by the application:

```bash
python -m pip install numpy pandas scipy scikit-learn matplotlib pyvista pyvistaqt PyQt5 pyserial openpyxl
```

Place the input files in the same folder as the script:

- `exportgeometry.csv` and `export1.csv` through `export42.csv` — ANSYS geometry and simulation exports.
- `DOexp.csv` — dissolved-oxygen model data.
- `inlet.stl`, `outlet.stl`, `wall.stl`, `symmetry.stl`, and `wall top.stl` — 3D boundary meshes.

### Run

```bash
python "Digital Twin Acq5.py"
```

The dashboard can be explored in test mode without an Arduino. Live readings use an Arduino connected on `COM3` at 9600 baud, sending temperature, flow, and pH as three comma-separated values per line.
change