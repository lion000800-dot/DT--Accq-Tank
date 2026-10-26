# Aquarium Digital Twin: Reduced-Order Model (ROM)

This repository contains a PyQt5 desktop application for exploring an aquarium digital twin. It reconstructs temperature and velocity fields from ANSYS simulation snapshots and estimates dissolved oxygen (DO) from the supplied DO data grid.

**A sensor is not required to run or test the ROM.** The application starts with default inputs and includes a slider-based test mode. A physical serial sensor is optional.

## What a reviewer needs

The application reads files by relative path from its current working directory. Start it from the project directory and make sure the following data files are beside `Digital Twin Acq4.py`.

> The simulation exports and meshes are not tracked in this source repository. A source-only clone cannot start the ROM until the authors provide these files separately. Please confirm that the reviewer has access to the files before asking them to run the application.

| Required file(s) | Format and use |
| --- | --- |
| `exportgeometry.csv` | ANSYS geometry export. It must contain `Node Number`, `X [ m ]`, `Y [ m ]`, and `Z [ m ]` columns. |
| `export1.csv` through `export42.csv` | Forty-two ANSYS snapshots. Each must contain `Node Number`, `Temperature [ C ]`, `Velocity u [ m s^-1 ]`, `Velocity v [ m s^-1 ]`, and `Velocity w [ m s^-1 ]`. Node numbers must match the geometry export. |
| `DOexp.csv` | Read without treating any row as a column header. The velocity axis is read from row 2, starting at column 3; the temperature axis is read from column 2, starting at row 3; grid values start at row 3, column 3. The axes and grid values must be numeric, ordered, and dimensionally consistent. Cubic interpolation requires at least four values on each axis. |
| `inlet.stl`, `outlet.stl`, `wall.stl`, `symmetry.stl`, `wall top.stl` | STL boundary meshes used to display the reconstructed fields. |

The application also looks for `icon_flow.png`, `icon_temp.png`, `icon_ph.png`, `icon_do.png`, `icon_warning.png`, and `SIITec logo.png` for the dashboard artwork. These are visual assets rather than ROM input data.

The other CSV and STL files in the project folder are not read by this application. Keep the column names, units, and node-number correspondence in the required files unchanged when preparing replacement exports.

## Install

The ROM engine was smoke-tested on Windows with Python 3.13.5. `requirements.txt` pins the direct dependency versions used for that check. Other Python and operating-system combinations have not been verified.

From the project directory, create an environment, install dependencies, and start the app:

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe "Digital Twin Acq4.py"
```

Run the final command from the project directory so the application can find the input files and artwork. The ROM setup reads all 42 snapshots and may take several seconds; the setup and prediction times are printed to the console.

## Run the ROM without a sensor

1. Start the app with the command above. Do not click the serial connection button. The app begins with default inputs of **25 °C**, **0 L/min**, and **pH 7** and calculates the initial fields without a sensor.
2. In the dashboard, choose **Temperature View**, **Velocity View**, or **DO View** to inspect the corresponding field. Choose **Probe View** to inspect the local temperature, velocity, and DO readouts.
3. Click **TEST** to show the test sliders. Set temperature to **25 °C** and flow to **2.5 L/min** (both are within the simulation training range), then move a slider. The reconstructed fields and probe values update automatically. The pH slider changes the displayed pH input; pH is not an input to the ROM field reconstruction.
4. Click **STOP TEST** to hide the sliders. The most recent test values remain active.

The temperature and flow snapshots span 24–30 °C and 0–5 L/min. The sliders extend beyond those ranges, but use inputs within the snapshot ranges when demonstrating interpolation. The DO estimate additionally depends on the local velocity/temperature prediction and the values in `DOexp.csv`.

The dashboard's test mode is an interactive demonstration, not an automated test suite or experimental validation. The validation window can generate points for analysis, but generated points are not experimental observations. Enter real experimental DO values in the **Exp DO (Input)** column before calculating errors; its default values are zero and do not provide meaningful validation metrics.

## Optional serial sensor

The app can read one comma-separated line containing numeric temperature, flow, and pH values:

```text
25.0,2.5,7.0
```

The serial connection is optional and uses 9600 baud. The app attempts to detect a serial device when the connection button is used. In test mode, clicking **TEST** disconnects an active serial connection.

## Reproduction and publication

For a reviewer-facing demonstration, provide the required simulation exports and meshes along with access to the source repository, then follow the sensor-free steps above. Do not present test-mode results as experimental evidence. See [JOURNAL_SUBMISSION.md](./JOURNAL_SUBMISSION.md) for the manuscript, software citation, data availability, and publication checklist.

This repository currently has no software license or automated test suite. Confirm permission to redistribute the simulation data and meshes, and add an approved license before describing the software as open source.
