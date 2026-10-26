# Digital Twin for Aquarium Monitoring

A desktop application for visualizing an aquarium digital twin. The application builds reduced-order models from ANSYS simulation exports and presents predicted temperature, velocity, and dissolved oxygen in an interactive 3D dashboard. It also includes a probe view, alarm indicators, a test mode, and an experimental dissolved-oxygen validation workflow with plots and Excel export.

## Repository contents

- `Digital Twin Acq5.py` — the PyQt5 dashboard application.
- `requirements.txt` — Python package dependencies.
- `JOURNAL_SUBMISSION.md` — a checklist for preparing a reproducible journal submission.

The simulation data and STL boundary meshes are not included in this repository. The application cannot start its full dashboard without these input files.

## Requirements

- Python 3 with the packages listed in `requirements.txt`.
- The input files listed below, in the same working directory from which the application is launched.
- Optional: an Arduino connected to `COM3` for live readings. The dashboard's test mode can be used without an Arduino.

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Input files

Place the following files beside the Python script:

| File(s) | Expected contents |
| --- | --- |
| `exportgeometry.csv` | ANSYS geometry table with `Node Number`, `X [ m ]`, `Y [ m ]`, and `Z [ m ]` columns. |
| `export1.csv` through `export42.csv` | ANSYS result tables, each with `Node Number`, `Temperature [ C ]`, `Velocity u [ m s^-1 ]`, `Velocity v [ m s^-1 ]`, and `Velocity w [ m s^-1 ]` columns. The node numbers must match the geometry file. |
| `DOexp.csv` | Numeric dissolved-oxygen data arranged as a grid: velocity values are read from row 2 starting at column 3; temperature values are read from column 2 starting at row 3; the corresponding measurements fill the grid from row 3, column 3. The file is read without a header row. |
| `inlet.stl`, `outlet.stl`, `wall.stl`, `symmetry.stl`, `wall top.stl` | STL boundary meshes used to display the model surfaces. |

The ANSYS loader locates the header row by searching for `Node` and then selects the named columns above. Use compatible units and data layouts when preparing exports.

## Run

From the directory containing the input files:

```bash
python "Digital Twin Acq5.py"
```

Use the dashboard's test mode to explore the application without connecting hardware. Live serial input is hard-coded to `COM3` at 9600 baud. Send one line per reading in this format:

```text
temperature,flow,pH
```

For example, `25.0,2.5,7.0`. Values must be numeric and comma-separated.

## Validation and export

The dissolved-oxygen validation window can generate sample points, accept experimental DO values, calculate error metrics, plot error distributions, and export the results to an Excel workbook. Experimental measurements are not included; provide them through the application before interpreting validation results.

## Reproducibility and publication

See [JOURNAL_SUBMISSION.md](./JOURNAL_SUBMISSION.md) for a checklist covering data availability, software citation, licensing, and manuscript reporting. Do not publish input data, meshes, or code unless you have permission to redistribute them.

## Limitations

- Input paths and the serial port are currently hard-coded; run the application from the directory containing its input files.
- Startup requires the ANSYS exports and STL meshes described above.
- This repository does not currently include a software license or a formal test suite. Confirm licensing and journal requirements before redistribution or publication.
