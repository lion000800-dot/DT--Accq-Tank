import sys
import time
import pandas as pd
import numpy as np
import pyvista as pv
import serial
from serial.tools import list_ports
import matplotlib.pyplot as plt 
from PyQt5 import QtWidgets, QtCore, QtGui
from pyvistaqt import QtInteractor
from scipy.linalg import svd
from scipy.interpolate import Rbf
from itertools import product
from scipy.spatial import KDTree
from scipy.interpolate import RegularGridInterpolator
from sklearn.linear_model import LinearRegression

# --- 0. DO MODEL SETUP ---
class DOModel:
    def __init__(self, csv_path='DOexp.csv'):
        try:
            df = pd.read_csv(csv_path, header=None)
            self.v_vals = df.iloc[1, 2:].astype(float).values    
            self.t_vals = df.iloc[2:, 1].astype(float).values 
            self.z_data = df.iloc[2:, 2:].astype(float).values   
            
            self.interpolator = RegularGridInterpolator(
                (self.t_vals, self.v_vals), 
                self.z_data, 
                method='cubic', 
                bounds_error=False, 
                fill_value=None  
            )
        except Exception as e:
            print(f"DO Model Load Error: {e}")
            self.interpolator = None 

    def predict_array(self, velocities, temperatures):
        if self.interpolator is None:
            return np.zeros_like(velocities)
        
        pts = np.column_stack((temperatures, velocities))
        results = self.interpolator(pts)
        return np.maximum(0, results) 
    
# --- 1. DATA LOADING UTILITY ---
def load_ansys_csv(path, cols):
    try:
        with open(path, 'r') as f: lines = f.readlines()
        start = next(i for i, line in enumerate(lines) if "Node" in line)
        df = pd.read_csv(path, skiprows=start)
        df.columns = [c.strip() for c in df.columns]
        return df[cols].apply(pd.to_numeric, errors='coerce').dropna()
    except: return None

# --- 2. DO VALIDATION WINDOW (UPDATED WITH LIVE RECALCULATION) ---
class DOValidationWindow(QtWidgets.QDialog):
    def __init__(self, main_app):
        super().__init__(main_app)
        self.main_app = main_app
        self.setWindowTitle("DO Surrogate Experimental Validation")
        self.resize(1000, 650)
        self.setStyleSheet("background-color: #f0f3f4;")
        
        self.plot_data = None 
        self.is_updating = False 
        
        layout = QtWidgets.QVBoxLayout(self)
        
        lbl_header = QtWidgets.QLabel("Random Sampling & Error Distribution")
        lbl_header.setStyleSheet("font-size: 20px; font-weight: bold; color: #2C3E50;")
        layout.addWidget(lbl_header)
        
        ctrl_layout = QtWidgets.QHBoxLayout()
        self.spin_n = QtWidgets.QSpinBox()
        self.spin_n.setRange(1, 100)
        self.spin_n.setValue(10)
        
        btn_generate = QtWidgets.QPushButton("1. Generate Points")
        btn_generate.setStyleSheet("background-color: #2980B9; color: white; padding: 8px; font-weight: bold;")
        btn_generate.clicked.connect(self.generate_points)
        
        btn_calc = QtWidgets.QPushButton("2. Calculate Errors")
        btn_calc.setStyleSheet("background-color: #27AE60; color: white; padding: 8px; font-weight: bold;")
        btn_calc.clicked.connect(self.calculate_errors)

        btn_plot = QtWidgets.QPushButton("3. Plot Error Distributions")
        btn_plot.setStyleSheet("background-color: #8E44AD; color: white; padding: 8px; font-weight: bold;")
        btn_plot.clicked.connect(self.show_plots)

        btn_export = QtWidgets.QPushButton("4. Export to Excel")
        btn_export.setStyleSheet("background-color: #E67E22; color: white; padding: 8px; font-weight: bold;")
        btn_export.clicked.connect(self.export_excel)
        
        ctrl_layout.addWidget(QtWidgets.QLabel("Random Points:"))
        ctrl_layout.addWidget(self.spin_n)
        ctrl_layout.addWidget(btn_generate)
        ctrl_layout.addWidget(btn_calc)
        ctrl_layout.addWidget(btn_plot)
        ctrl_layout.addWidget(btn_export)
        ctrl_layout.addStretch()
        layout.addLayout(ctrl_layout)
        
        self.table = QtWidgets.QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(["Q (L/min)", "T (°C)", "X (m)", "Y (m)", "Z (m)", "ROM DO (Pred)", "Exp DO (Input)", "Local Error %"])
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.table.setStyleSheet("background-color: white; gridline-color: #d0d3d4;")
        
        self.table.cellChanged.connect(self.on_cell_changed)
        layout.addWidget(self.table)
        
        res_layout = QtWidgets.QHBoxLayout()
        self.lbl_mae = QtWidgets.QLabel("MAE: -- mg/L")
        self.lbl_rmse = QtWidgets.QLabel("RMSE: -- mg/L")
        self.lbl_global = QtWidgets.QLabel("Global Relative Error: -- %")
        
        for lbl in [self.lbl_mae, self.lbl_rmse, self.lbl_global]:
            lbl.setStyleSheet("font-size: 16px; font-weight: bold; color: #C0392B; padding: 10px; border: 1px solid #C0392B; background-color: #FDEDEC;")
            res_layout.addWidget(lbl)
            
        layout.addLayout(res_layout)

    def generate_points(self):
        self.is_updating = True 

        n = self.spin_n.value()
        self.table.setRowCount(n)
        
        Q_rand = np.round(np.random.uniform(2.0, 5.0, n), 1)
        T_rand = np.round(np.random.uniform(24.0, 30.0, n), 1)
        X_rand = np.round(np.random.uniform(0.0, 0.39, n), 2)
        Y_rand = np.round(np.random.uniform(0.0, 0.15, n), 2)
        Z_rand = np.round(np.random.uniform(-0.10, 0.0, n), 2)
        
        for i in range(n):
            q, t = Q_rand[i], T_rand[i]
            x, y, z = X_rand[i], Y_rand[i], Z_rand[i]
            
            a_T = np.array([f(t, q) for f in self.main_app.rbf_T])
            a_U = np.array([f(t, q) for f in self.main_app.rbf_U])
            
            T_recon = self.main_app.Phi_T @ a_T
            U_recon = (self.main_app.Phi_U @ a_U).reshape(3, self.main_app.N).T
            Vmag = np.linalg.norm(U_recon, axis=1)
            
            _, idx = self.main_app.rom_tree.query([x, y, z])
            T_p, V_p = T_recon[idx], Vmag[idx]
            DO_pred = self.main_app.solver_do.predict_array(np.array([V_p]), np.array([T_p]))[0]
            
            self.table.setItem(i, 0, QtWidgets.QTableWidgetItem(f"{q:.1f}"))
            self.table.setItem(i, 1, QtWidgets.QTableWidgetItem(f"{t:.1f}"))
            self.table.setItem(i, 2, QtWidgets.QTableWidgetItem(f"{x:.2f}"))
            self.table.setItem(i, 3, QtWidgets.QTableWidgetItem(f"{y:.2f}"))
            self.table.setItem(i, 4, QtWidgets.QTableWidgetItem(f"{z:.2f}"))
            
            pred_item = QtWidgets.QTableWidgetItem(f"{DO_pred:.3f}")
            pred_item.setFlags(QtCore.Qt.ItemIsEnabled) 
            self.table.setItem(i, 5, pred_item)
            
            exp_item = QtWidgets.QTableWidgetItem("0.00")
            exp_item.setBackground(QtGui.QColor("#FEF9E7")) 
            self.table.setItem(i, 6, exp_item)
            
            err_item = QtWidgets.QTableWidgetItem("--")
            err_item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(i, 7, err_item)

        self.plot_data = None 
        self.is_updating = False 

    def on_cell_changed(self, row, col):
        if self.is_updating:
            return 

        if col in [0, 1, 2, 3, 4]:
            try:
                q = float(self.table.item(row, 0).text())
                t = float(self.table.item(row, 1).text())
                x = float(self.table.item(row, 2).text())
                y = float(self.table.item(row, 3).text())
                z = float(self.table.item(row, 4).text())
                
                a_T = np.array([f(t, q) for f in self.main_app.rbf_T])
                a_U = np.array([f(t, q) for f in self.main_app.rbf_U])
                
                T_recon = self.main_app.Phi_T @ a_T
                U_recon = (self.main_app.Phi_U @ a_U).reshape(3, self.main_app.N).T
                Vmag = np.linalg.norm(U_recon, axis=1)
                
                _, idx = self.main_app.rom_tree.query([x, y, z])
                T_p, V_p = T_recon[idx], Vmag[idx]
                DO_pred = self.main_app.solver_do.predict_array(np.array([V_p]), np.array([T_p]))[0]
                
                self.is_updating = True
                self.table.item(row, 5).setText(f"{DO_pred:.3f}")
                
                if self.table.item(row, 7) is not None:
                    self.table.item(row, 7).setText("--")
                    
                self.is_updating = False
                
            except (ValueError, AttributeError):
                pass

    def calculate_errors(self):
        n = self.table.rowCount()
        data = []
        
        for i in range(n):
            try:
                q = float(self.table.item(i, 0).text())
                t = float(self.table.item(i, 1).text())
                x = float(self.table.item(i, 2).text())
                y = float(self.table.item(i, 3).text())
                z = float(self.table.item(i, 4).text())
                pred = float(self.table.item(i, 5).text())
                exp = float(self.table.item(i, 6).text())
                
                if exp > 0:
                    err = abs(exp - pred) / exp * 100
                else:
                    err = 0.0
                
                self.is_updating = True 
                self.table.setItem(i, 7, QtWidgets.QTableWidgetItem(f"{err:.2f}"))
                self.is_updating = False
                
                data.append({'Point ID': i+1, 'Q (L/min)': q, 'T (°C)': t, 'X (m)': x, 'Y (m)': y, 'Z (m)': z, 
                             'ROM DO': pred, 'Exp DO': exp, 'Local Error (%)': err})
            except ValueError:
                QtWidgets.QMessageBox.warning(self, "Input Error", f"Invalid input at row {i+1}")
                return

        self.plot_data = pd.DataFrame(data)
        
        preds = self.plot_data['ROM DO'].values
        exps = self.plot_data['Exp DO'].values
        
        mae = np.mean(np.abs(exps - preds))
        rmse = np.sqrt(np.mean((exps - preds)**2))
        
        norm_exp = np.linalg.norm(exps)
        global_err = (np.linalg.norm(exps - preds) / norm_exp * 100) if norm_exp > 0 else 0.0
            
        self.lbl_mae.setText(f"MAE: {mae:.4f} mg/L")
        self.lbl_rmse.setText(f"RMSE: {rmse:.4f} mg/L")
        self.lbl_global.setText(f"Global Relative Error: {global_err:.2f} %")

    def show_plots(self):
        if self.plot_data is None:
            QtWidgets.QMessageBox.warning(self, "Warning", "Please calculate errors first!")
            return
            
        df = self.plot_data
        
        fig, axs = plt.subplots(1, 3, figsize=(16, 5))
        fig.canvas.manager.set_window_title('DO Surrogate Error Analysis')
        
        axs[0].bar(df['Point ID'], df['Local Error (%)'], color='skyblue', edgecolor='black')
        axs[0].set_title('Local Relative Error per Spatial Point', fontweight='bold')
        axs[0].set_xlabel('Random Test Point ID')
        axs[0].set_ylabel('Relative Error (%)')
        axs[0].set_xticks(df['Point ID'])
        
        axs[1].scatter(df['Q (L/min)'], df['Local Error (%)'], color='coral', s=100, edgecolors='black', alpha=0.8)
        axs[1].set_title('Error Distribution vs. Flow Rate', fontweight='bold')
        axs[1].set_xlabel('Flow Rate (L/min)')
        axs[1].set_ylabel('Relative Error (%)')
        axs[1].grid(True, linestyle='--', alpha=0.7)
        
        axs[2].scatter(df['T (°C)'], df['Local Error (%)'], color='mediumseagreen', s=100, edgecolors='black', alpha=0.8)
        axs[2].set_title('Error Distribution vs. Temperature', fontweight='bold')
        axs[2].set_xlabel('Temperature (°C)')
        axs[2].set_ylabel('Relative Error (%)')
        axs[2].grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.show()

    def export_excel(self):
        if self.plot_data is None:
            QtWidgets.QMessageBox.warning(self, "Warning", "Please calculate errors first!")
            return
            
        options = QtWidgets.QFileDialog.Options()
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Validation Results", "DO_Validation_Results.xlsx", "Excel Files (*.xlsx)", options=options)
        
        if file_path:
            try:
                metrics = pd.DataFrame([{
                    'Point ID': 'METRICS',
                    'ROM DO': self.lbl_mae.text(),
                    'Exp DO': self.lbl_rmse.text(),
                    'Local Error (%)': self.lbl_global.text()
                }])
                df_export = pd.concat([self.plot_data, metrics], ignore_index=True)
                
                df_export.to_excel(file_path, index=False)
                QtWidgets.QMessageBox.information(self, "Success", f"Validation data successfully exported to:\n{file_path}")
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Export Error", f"Failed to save Excel file:\n{str(e)}")

# --- 3. MAIN APPLICATION ---
class ROMApp(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Digital Twin Acq - Live High-Visibility Dashboard")
        self.display_mode = "Temperature" 
        self.camera_view = "XY"
        self.vector_submode = "Top"
        self.pt = np.array([0.0, 0.0, 0.0])
        self.ser = None
        self.connected_port = None
        self.current_inputs = {'temp': 25.0, 'flow': 0.0, 'ph': 7.0}
        
        self.is_test_mode = False 

        # ตัวแปรสำหรับระบบแจ้งเตือน (Alarm & Blinking)
        self.alarms = {'temp': False, 'flow': False, 'ph': False, 'do': False}
        self.blink_state = False
        self.blink_timer = QtCore.QTimer()
        self.blink_timer.timeout.connect(self.blink_warnings)
        self.blink_timer.start(500) # กระพริบทุกๆ 0.5 วินาที

        self.setup_engine()
        self.setup_ui()
        
        self.update_connection_ui_state(False)

        self.acq_timer = QtCore.QTimer()
        self.acq_timer.timeout.connect(self.read_serial_data)
        self.acq_timer.start(500) 

    def setup_engine(self):
        start_time = time.perf_counter()
        T_axis = np.linspace(24, 30, 7); Q_axis = np.linspace(0, 5, 6)
        self.params_train = np.array(list(product(T_axis, Q_axis))) 
        df_geo = load_ansys_csv('exportgeometry.csv', ['Node Number', 'X [ m ]', 'Y [ m ]', 'Z [ m ]'])
        self.snaps = {'T': [], 'u': [], 'v': [], 'w': []}
        for i in range(1, 43):
            df_p = load_ansys_csv(f'export{i}.csv', ['Node Number', 'Temperature [ C ]', 'Velocity u [ m s^-1 ]', 'Velocity v [ m s^-1 ]', 'Velocity w [ m s^-1 ]'])
            m = pd.merge(df_geo, df_p, on='Node Number').sort_values('Node Number')
            for k in self.snaps: self.snaps[k].append(m['Temperature [ C ]' if k=='T' else f'Velocity {k} [ m s^-1 ]'].values)
        self.r = 8; self.N = len(self.snaps['T'][0])
        S_T = np.array(self.snaps['T']).T; S_U = np.vstack([np.array(self.snaps['u']).T, np.array(self.snaps['v']).T, np.array(self.snaps['w']).T])
        U_T, _, _ = svd(S_T, full_matrices=False); self.Phi_T = U_T[:, :self.r]; Alpha_T = self.Phi_T.T @ S_T
        self.rbf_T = [Rbf(self.params_train[:, 0], self.params_train[:, 1], Alpha_T[j, :], function='multiquadric') for j in range(self.r)]
        U_U, _, _ = svd(S_U, full_matrices=False); self.Phi_U = U_U[:, :self.r]; Alpha_U = self.Phi_U.T @ S_U
        self.rbf_U = [Rbf(self.params_train[:, 0], self.params_train[:, 1], Alpha_U[j, :], function='multiquadric') for j in range(self.r)]
        self.solver_do = DOModel('DOexp.csv'); self.coords = m[['X [ m ]', 'Y [ m ]', 'Z [ m ]']].values
        self.cloud = pv.PolyData(self.coords); self.rom_tree = KDTree(self.coords)
        self.boundaries = {}
        for name in ['inlet', 'outlet', 'wall', 'symmetry', 'wall top']:
            mesh = pv.read(f"{name}.stl")
            _, idx = self.rom_tree.query(mesh.points); self.boundaries[name] = {'mesh': mesh, 'indices': idx}
        self.clim_temp = [S_T.min(), S_T.max()]; self.clim_vel = [0, np.sqrt(S_U**2).max()]; self.clim_Do = [3, 10]
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"ROM setup time: {elapsed_time:.4f} seconds")

    def setup_ui(self):
        cw = QtWidgets.QWidget()
        self.setCentralWidget(cw)
        cw.setStyleSheet("background-color: #e6e6e6;") 
        
        main_layout = QtWidgets.QVBoxLayout(cw)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(8) 

        # --- 1. HEADER ---
        header = QtWidgets.QLabel("Digital Twin for Aquarium Monitoring")
        header.setStyleSheet("background-color: #1A5276; color: white; font-size: 35px; font-weight: bold; padding: 15px;")
        header.setAlignment(QtCore.Qt.AlignCenter)
        main_layout.addWidget(header)

        # Content Layout
        content_layout = QtWidgets.QHBoxLayout()
        content_layout.setContentsMargins(35, 10, 35, 20)
        main_layout.addLayout(content_layout)

        # ==========================================
        #               ฝั่งซ้าย (3D & Controls)
        # ==========================================
        left_layout = QtWidgets.QVBoxLayout()
        content_layout.addLayout(left_layout, stretch=7)

        # 1.1 Flow Rate Banner
        flow_frame = QtWidgets.QFrame()
        flow_frame.setStyleSheet("background-color: #1A5276; border-radius: 0px;")
        flow_layout = QtWidgets.QHBoxLayout(flow_frame)
        flow_layout.setContentsMargins(10, 8, 10, 8)

        icon_flow_lbl = QtWidgets.QLabel()
        pixmap_flow = QtGui.QPixmap("icon_flow.png")
        pixmap_flow = pixmap_flow.scaled(100, 100, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        icon_flow_lbl.setPixmap(pixmap_flow)
        
        flow_title = QtWidgets.QLabel(" Flow Rate [L/min]:")
        flow_title.setStyleSheet("color: white; font-size: 25px; font-weight: bold; border: none;")
        
        white_box_flow = QtWidgets.QFrame()
        white_box_flow.setStyleSheet("background-color: white; border: 1px solid black;")
        white_box_flow.setFixedWidth(150)  
        
        wb_flow_layout = QtWidgets.QHBoxLayout(white_box_flow)
        wb_flow_layout.setContentsMargins(5, 2, 10, 2) 
        wb_flow_layout.setSpacing(5)

        self.card_q_in = QtWidgets.QLabel("0.00")
        self.card_q_in.setStyleSheet("color: black; font-size: 25px; font-weight: bold; border: none;")
        self.card_q_in.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)

        self.warn_q_in = QtWidgets.QLabel()
        warn_pix_q = QtGui.QPixmap("icon_warning.png").scaled(25, 25, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        self.warn_q_in.setPixmap(warn_pix_q)
        self.warn_q_in.setFixedSize(25, 25)
        self.warn_q_in.setVisible(False) 
        self.warn_q_in.setStyleSheet("border: none;")

        wb_flow_layout.addWidget(self.warn_q_in)
        wb_flow_layout.addWidget(self.card_q_in, stretch=1)
        
        flow_layout.addStretch()
        flow_layout.addWidget(icon_flow_lbl) 
        flow_layout.addWidget(flow_title)
        flow_layout.addWidget(white_box_flow) 
        flow_layout.addStretch()
        left_layout.addWidget(flow_frame)

        # 1.2 3D Plotter
        self.plotter = QtInteractor(self, shape=(1, 1))
        self.plotter.set_background('white') 
        left_layout.addWidget(self.plotter.interactor, stretch=1)

        # ==========================================
        # 1.3 กล่อง TEST MODE (Sliders)
        # ==========================================
        self.test_panel = QtWidgets.QGroupBox("TEST MODE CONTROLS")
        self.test_panel.setStyleSheet("QGroupBox { border: 2px dashed #3498DB; margin-top: 10px; font-size: 15px; background-color: white;} QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px; color: #3498DB; font-weight: bold;}")
        test_vbox = QtWidgets.QVBoxLayout()
        test_vbox.setContentsMargins(15, 15, 15, 15)
        
        def create_test_slider(layout, label_text, min_val, max_val, start_val, param_key):
            h_layout = QtWidgets.QHBoxLayout()
            lbl = QtWidgets.QLabel(label_text)
            lbl.setFixedWidth(120)
            lbl.setStyleSheet("font-size: 14px; font-weight: bold;")

            val_lbl = QtWidgets.QLabel(f"{start_val:.2f}")
            val_lbl.setFixedWidth(50)
            val_lbl.setStyleSheet("font-size: 14px; font-weight: bold; color: #2980B9;")

            slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            slider.setMinimum(int(min_val * 100))
            slider.setMaximum(int(max_val * 100))
            slider.setValue(int(start_val * 100))

            def on_change(val):
                real_val = val / 100.0
                val_lbl.setText(f"{real_val:.2f}")
                self.test_val_changed(param_key, real_val)

            slider.valueChanged.connect(on_change)
            h_layout.addWidget(lbl)
            h_layout.addWidget(slider)
            h_layout.addWidget(val_lbl)
            layout.addLayout(h_layout)
            return slider

        self.slider_temp = create_test_slider(test_vbox, "Temp [°C]:", 23.0, 31.0, 25.0, 'temp')
        self.slider_flow = create_test_slider(test_vbox, "Flow [L/min]:", 0.0, 6.0, 0.0, 'flow')
        self.slider_ph = create_test_slider(test_vbox, "pH:", 0.0, 14.0, 7.0, 'ph')

        self.test_panel.setLayout(test_vbox)
        self.test_panel.setVisible(False) 
        left_layout.addWidget(self.test_panel)

        # 1.4 Bottom Controls
        bottom_ctrl_layout = QtWidgets.QHBoxLayout()
        
        group_display = QtWidgets.QGroupBox("DISPLAY VIEW")
        group_display.setStyleSheet("QGroupBox { border: 1px solid gray; margin-top: 10px; font-size: 15px; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }")
        
        disp_grid = QtWidgets.QGridLayout()
        views = [("Velocity View", "Velocity"), ("Temperature View", "Temperature"), ("DO View", "Do"), ("Probe View", "Probe"), ("Vectors View", "Vectors")]
        
        for i, (txt, mode) in enumerate(views):
            btn = QtWidgets.QPushButton(txt)
            btn.setStyleSheet("background-color: white; border: 1px solid gray; padding: 5px; font-size: 14px;")
            btn.clicked.connect(lambda checked, m=mode: self.set_view(m))
            
            if i == 4: disp_grid.addWidget(btn, 2, 0, 1, 2) 
            else: disp_grid.addWidget(btn, i//2, i%2)
                
        group_display.setLayout(disp_grid)
        bottom_ctrl_layout.addWidget(group_display)

        self.vector_panel = QtWidgets.QGroupBox("Vector Orientation")
        self.vector_panel.setStyleSheet("QGroupBox { border: 1px solid gray; margin-top: 10px; font-size: 15px; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }")

        vec_vbox = QtWidgets.QVBoxLayout()
        btn_top = QtWidgets.QPushButton("TOP (XY)")
        btn_top.setStyleSheet("background-color: white; border: 1px solid gray; padding: 5px; font-size: 14px;")
        btn_top.clicked.connect(lambda: self.set_vector_view("Top"))
        btn_side = QtWidgets.QPushButton("Side (XY)")
        btn_side.setStyleSheet("background-color: white; border: 1px solid gray; padding: 5px; font-size: 14px;")
        btn_side.clicked.connect(lambda: self.set_vector_view("Side"))
        vec_vbox.addWidget(btn_top); vec_vbox.addWidget(btn_side)
        self.vector_panel.setLayout(vec_vbox)
        self.vector_panel.setVisible(False) 
        bottom_ctrl_layout.addWidget(self.vector_panel)

        group_cam = QtWidgets.QGroupBox("CAMERA ORIENTATION")
        group_cam.setStyleSheet("QGroupBox { border: 1px solid gray; margin-top: 10px; font-size: 15px; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }")

        cam_grid = QtWidgets.QGridLayout()
        cams = [("Front (XY)", "XY"), ("Isometric view", "ISO"), ("Top (XZ)", "XZ")]
        for i, (txt, mode) in enumerate(cams):
            btn = QtWidgets.QPushButton(txt)
            btn.setStyleSheet("background-color: white; border: 1px solid gray; padding: 5px; font-size: 14px;")
            btn.clicked.connect(lambda checked, m=mode: self.set_camera(m))
            cam_grid.addWidget(btn, i//2, i%2)
        group_cam.setLayout(cam_grid)
        bottom_ctrl_layout.addWidget(group_cam)

        left_layout.addLayout(bottom_ctrl_layout)

        # ==========================================
        #               ฝั่งขวา (Sensors & Probe)
        # ==========================================
        right_layout = QtWidgets.QVBoxLayout()
        content_layout.addLayout(right_layout, stretch=2) 

        # 2.1 Sensor Inputs
        sensor_box = QtWidgets.QVBoxLayout()
        self.card_t_in, self.warn_t_in = self.create_right_sensor(sensor_box, "  Temp [°C ]:", "icon_temp.png")
        self.card_ph, self.warn_ph = self.create_right_sensor(sensor_box, "  pH:", "icon_ph.png")
        self.card_do_in, self.warn_do = self.create_right_sensor(sensor_box, " DO [mg/L]:", "icon_do.png")
        
        right_layout.addLayout(sensor_box)

        # 2.2 Logo Placeholder (ใส่โลโก้คณะแทนกล่องว่าง)
        self.graph_frame = QtWidgets.QFrame()
        self.graph_frame.setStyleSheet("background-color: white; border: 1px solid #CCC;")
        self.graph_frame.setMinimumHeight(200)
        
        # สร้าง Layout ภายในกล่องสีขาว
        logo_layout = QtWidgets.QVBoxLayout(self.graph_frame)
        logo_layout.setContentsMargins(10, 10, 10, 10)
        
        # โหลดรูปภาพโลโก้
        self.logo_lbl = QtWidgets.QLabel()
        pixmap_logo = QtGui.QPixmap("SIITec logo.png") 
        
        # ปรับขนาดรูปให้พอดี (สามารถเปลี่ยนเลข 300 เพื่อปรับความเล็ก/ใหญ่ได้ตามต้องการ)
        pixmap_logo = pixmap_logo.scaled(300, 300, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        
        self.logo_lbl.setPixmap(pixmap_logo)
        self.logo_lbl.setAlignment(QtCore.Qt.AlignCenter) # จัดรูปให้อยู่กึ่งกลางกล่อง
        
        logo_layout.addWidget(self.logo_lbl)
        right_layout.addWidget(self.graph_frame, stretch=1)

        # 2.3 Probe Position Group
        group_probe = QtWidgets.QGroupBox("PROBE POSITION [m]")
        group_probe.setStyleSheet("QGroupBox { border: 1px solid gray; margin-top: 10px; font-size: 15px; } QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }")
        probe_vbox = QtWidgets.QVBoxLayout()
        
        self.in_x = self.create_input_field(probe_vbox, "  X:")
        self.in_y = self.create_input_field(probe_vbox, "  Y:")
        self.in_z = self.create_input_field(probe_vbox, "  Z:")
        
        btn_probe_layout = QtWidgets.QHBoxLayout()
        btn_probe_update = QtWidgets.QPushButton("UPDATE PROBE POSITION")
        btn_probe_update.setStyleSheet("background-color: #7D3C98; color: white; font-weight: bold; padding: 10px; margin: 10px 0px;")
        btn_probe_update.clicked.connect(self.manual_probe_sync)
        
        self.btn_test = QtWidgets.QPushButton("TEST")
        self.btn_test.setStyleSheet("background-color: #3498DB; color: white; font-weight: bold; padding: 10px; margin: 10px 0px;")
        self.btn_test.clicked.connect(self.toggle_test_mode)
        
        btn_probe_layout.addWidget(btn_probe_update, stretch=7) 
        btn_probe_layout.addWidget(self.btn_test, stretch=3)    
        probe_vbox.addLayout(btn_probe_layout)

        self.card_res_v = self.create_output_field(probe_vbox, "Probe Velocity [m/s]:")
        self.card_res_t = self.create_output_field(probe_vbox, "Probe Temp [°C]:")
        self.card_res_do = self.create_output_field(probe_vbox, "Probe DO [mg/L]:")
        
        group_probe.setLayout(probe_vbox)
        right_layout.addWidget(group_probe)

        # --- เพิ่มปุ่ม OPEN DO VALIDATION MODULE ตรงนี้ตามรูปเป๊ะๆ ---
        self.btn_open_val = QtWidgets.QPushButton("OPEN DO VALIDATION MODULE")
        self.btn_open_val.setStyleSheet("background-color: #E67E22; color: white; font-weight: bold; padding: 10px; margin: 5px 0px; border: 1px solid black;")
        self.btn_open_val.clicked.connect(self.open_validation_window)
        right_layout.addWidget(self.btn_open_val)
        # --------------------------------------------------------

        # 2.4 Connect Buttons
        connection_layout = QtWidgets.QHBoxLayout()
        self.btn_connect_com = QtWidgets.QPushButton("CONNECTED (COM)")
        self.btn_connect_com.clicked.connect(self.connect_serial)
        self.btn_disconnect = QtWidgets.QPushButton("DISCONNECTED")
        self.btn_disconnect.clicked.connect(self.disconnect_serial)
        
        connection_layout.addWidget(self.btn_connect_com)
        connection_layout.addWidget(self.btn_disconnect)
        right_layout.addLayout(connection_layout)

        self.update_plots()

    # --- UI Helper Functions ---
    def create_right_sensor(self, layout, title, icon_path=None):
        frame = QtWidgets.QFrame()
        frame.setStyleSheet("background-color: #1A5276; border-radius: 0px;")
        h_layout = QtWidgets.QHBoxLayout(frame)
        h_layout.setContentsMargins(10, 8, 10, 8)
        
        if icon_path:
            icon_lbl = QtWidgets.QLabel()
            pixmap = QtGui.QPixmap(icon_path)
            pixmap = pixmap.scaled(37, 37, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
            icon_lbl.setPixmap(pixmap)
            h_layout.addWidget(icon_lbl)

        lbl = QtWidgets.QLabel(title)
        lbl.setStyleSheet("color: white; font-size: 25px; font-weight: bold; border: none;")
        
        # --- ปรับโครงสร้างกล่องสีขาว ฝั่งขวา (ย้ายไอคอนเตือนมาฝั่งซ้าย) ---
        white_box = QtWidgets.QFrame()
        white_box.setStyleSheet("background-color: white; border: 1px solid black;")
        white_box.setMinimumWidth(160)
        wb_layout = QtWidgets.QHBoxLayout(white_box)
        # ปรับขอบ: ซ้าย 5, บน 2, ขวา 10, ล่าง 2
        wb_layout.setContentsMargins(5, 2, 10, 2) 
        wb_layout.setSpacing(5)

        val_lbl = QtWidgets.QLabel("0.00")
        val_lbl.setStyleSheet("color: black; font-size: 25px; font-weight: bold; border: none;")
        val_lbl.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)

        warn_lbl = QtWidgets.QLabel()
        warn_pix = QtGui.QPixmap("icon_warning.png").scaled(25, 25, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        warn_lbl.setPixmap(warn_pix)
        warn_lbl.setFixedSize(25, 25)
        warn_lbl.setVisible(False) 
        warn_lbl.setStyleSheet("border: none;")

        # นำไอคอนเตือนเข้าก่อน เพื่อให้อยู่ซ้ายมือ แล้วค่อยตามด้วยตัวเลข
        wb_layout.addWidget(warn_lbl)
        wb_layout.addWidget(val_lbl, stretch=1)
        # --------------------------------------

        h_layout.addWidget(lbl)
        h_layout.addStretch()
        h_layout.addWidget(white_box)
        layout.addWidget(frame)
        
        return val_lbl, warn_lbl

    def create_input_field(self, layout, label_text):
        h_layout = QtWidgets.QHBoxLayout()
        lbl = QtWidgets.QLabel(label_text)
        lbl.setStyleSheet("font-size: 15px; font-weight: normal; color: #000000;")
        lbl.setFixedWidth(30)
        line_edit = QtWidgets.QLineEdit("0.00")
        line_edit.setStyleSheet("background-color: white; border: 1px solid gray; padding: 3px; font-size: 14px;")
        h_layout.addWidget(lbl)
        h_layout.addWidget(line_edit)
        layout.addLayout(h_layout)
        return line_edit

    def create_output_field(self, layout, label_text):
        h_layout = QtWidgets.QHBoxLayout()
        lbl = QtWidgets.QLabel(label_text)
        lbl.setStyleSheet("font-size: 15px; font-weight: normal; color: #000000;")

        lbl.setFixedWidth(150)
        val_lbl = QtWidgets.QLineEdit("0.00") 
        val_lbl.setReadOnly(True)
        val_lbl.setStyleSheet("background-color: white; border: 1px solid gray; padding: 3px; font-size: 14px;")
        h_layout.addWidget(lbl)
        h_layout.addWidget(val_lbl)
        layout.addLayout(h_layout)
        return val_lbl

    # ==========================================
    # ระบบการแจ้งเตือน (Alarm Systems)
    # ==========================================
    def update_alarms(self, t, q, ph, do):
        is_active = (self.ser is not None) or self.is_test_mode
        
        if is_active:
            self.alarms['flow'] = not (2.0 <= q <= 5.0)
            self.alarms['temp'] = not (24.0 <= t <= 30.0) 
            self.alarms['ph'] = not (7.0 <= ph <= 8.0)
            self.alarms['do'] = not (do >= 5.0)
        else:
            self.alarms = {k: False for k in self.alarms}

        # ฟังก์ชันช่วยเปลี่ยนสีตัวอักษร และจัดการซ่อนไอคอน
        def set_alarm_ui(lbl, warn_icon, is_alarm):
            color = "red" if is_alarm else "black"
            lbl.setStyleSheet(f"color: {color}; font-size: 25px; font-weight: bold; border: none;")
            
            # --- แก้ไขจุดที่ไอคอนค้าง ---
            if not is_alarm:
                warn_icon.setVisible(False)

        # สั่งอัปเดตทั้งตัวเลขและไอคอนไปพร้อมๆ กัน
        set_alarm_ui(self.card_q_in, self.warn_q_in, self.alarms['flow'])
        set_alarm_ui(self.card_t_in, self.warn_t_in, self.alarms['temp'])
        set_alarm_ui(self.card_ph, self.warn_ph, self.alarms['ph'])
        set_alarm_ui(self.card_do_in, self.warn_do, self.alarms['do'])

    def blink_warnings(self):
        self.blink_state = not self.blink_state
        
        # กระพริบเฉพาะค่าที่กำลังมีปัญหา (Alarm = True)
        if self.alarms['flow']: self.warn_q_in.setVisible(self.blink_state)
        if self.alarms['temp']: self.warn_t_in.setVisible(self.blink_state)
        if self.alarms['ph']: self.warn_ph.setVisible(self.blink_state)
        if self.alarms['do']: self.warn_do.setVisible(self.blink_state)
    # ==========================================

    # --- Action Functions ---
    def open_validation_window(self):
        """เปิดหน้าต่าง DO Validation Module"""
        val_win = DOValidationWindow(self)
        val_win.exec_()
        
    def toggle_test_mode(self):
        self.is_test_mode = not self.is_test_mode
        if self.is_test_mode:
            if self.ser is not None: self.disconnect_serial()
            self.test_panel.setVisible(True)
            self.btn_test.setText("STOP TEST")
            self.btn_test.setStyleSheet("background-color: #E74C3C; color: white; font-weight: bold; padding: 10px; margin: 10px 0px;") 
            self.test_val_changed('temp', self.slider_temp.value() / 100.0)
            self.test_val_changed('flow', self.slider_flow.value() / 100.0)
            self.test_val_changed('ph', self.slider_ph.value() / 100.0)
        else:
            self.test_panel.setVisible(False)
            self.btn_test.setText("TEST")
            self.btn_test.setStyleSheet("background-color: #3498DB; color: white; font-weight: bold; padding: 10px; margin: 10px 0px;") 
            self.update_plots()

    def test_val_changed(self, param, value):
        if self.is_test_mode:
            self.current_inputs[param] = value
            self.card_t_in.setText(f"{self.current_inputs['temp']:.2f}")
            self.card_q_in.setText(f"{self.current_inputs['flow']:.2f}")
            self.card_ph.setText(f"{self.current_inputs['ph']:.2f}")
            self.update_plots()

    def detect_serial_port(self):
        ports = list(list_ports.comports())
        if not ports:
            return None

        port_names = [p.device for p in ports]
        candidates = []
        for p in ports:
            desc = (p.description or "").lower()
            if any(key in desc for key in ["arduino", "usb serial", "cp210", "ch340", "ftdi", "wchusbserial"]):
                candidates.append(p.device)

        return candidates[0] if candidates else port_names[0]

    def connect_serial(self):
        if self.is_test_mode: self.toggle_test_mode()
        if self.ser is None:
            try:
                port_name = self.detect_serial_port()
                if port_name is None:
                    QtWidgets.QMessageBox.warning(self, "Connection Error", "No serial device detected.")
                    self.update_connection_ui_state(False)
                    return

                self.ser = serial.Serial(port_name, 9600, timeout=0.1)
                self.connected_port = port_name
                self.btn_connect_com.setText(f"CONNECTED ({port_name})")
                self.update_connection_ui_state(True)
            except Exception as exc:
                self.connected_port = None
                QtWidgets.QMessageBox.warning(self, "Connection Error", f"Arduino not detected.\n{exc}")
                self.update_connection_ui_state(False)

    def disconnect_serial(self):
        if self.ser is not None:
            self.ser.close()
            self.ser = None
        self.btn_connect_com.setText("CONNECTED (COM)")
        self.update_connection_ui_state(False)
        self.update_plots() 

    def update_connection_ui_state(self, is_connected):
        base_style = "color: white; font-weight: bold; padding: 10px; font-size: 14px;"
        if is_connected:
            self.btn_connect_com.setStyleSheet(f"background-color: #27AE60; {base_style}") 
            self.btn_disconnect.setStyleSheet(f"background-color: #95A5A6; {base_style}")  
            self.btn_connect_com.setEnabled(False) 
            self.btn_disconnect.setEnabled(True)   
        else:
            self.btn_connect_com.setStyleSheet(f"background-color: #95A5A6; {base_style}") 
            self.btn_disconnect.setStyleSheet(f"background-color: #E74C3C; {base_style}")  
            self.btn_connect_com.setEnabled(True)  
            self.btn_disconnect.setEnabled(False) 

    def read_serial_data(self):
        if self.ser and self.ser.in_waiting > 0:
            try:
                line = self.ser.readline().decode('utf-8', errors='ignore').strip()
                parts = line.split(',')
                if len(parts) == 3:
                    self.current_inputs = {'temp': float(parts[0]), 'flow': float(parts[1]), 'ph': float(parts[2])}
                    self.card_t_in.setText(f"{self.current_inputs['temp']:.1f}")
                    self.card_q_in.setText(f"{self.current_inputs['flow']:.2f}")
                    self.card_ph.setText(f"{self.current_inputs['ph']:.2f}")
                    self.update_plots() 
            except: pass

    def set_camera(self, mode):
        if mode == "XY": 
            self.plotter.view_xy()
        elif mode == "XZ": 
            self.plotter.view_xz(negative=True)
            self.plotter.camera.roll += 180
        else: 
            # ==========================================
            # --- สร้างมุมมอง Custom Isometric ด้วยตัวเอง ---
            # ==========================================
            self.plotter.view_xy()  # 1. เริ่มตั้งกล้องมองจากด้านหน้าตรงๆ
            
            # 2. หมุนกล้องแนวนอน (Azimuth) ให้เห็นด้านข้าง
            # (ถ้าหมุนแล้วสลับฝั่ง ให้ลองเปลี่ยนจาก += เป็น -= ดูนะครับ)
            self.plotter.camera.azimuth -= 25   
            
            # 3. ยกกล้องขึ้นแนวตั้ง (Elevation) ให้เห็นด้านบนนิดหน่อย
            self.plotter.camera.elevation += 15 
            # ==========================================
            
        self.plotter.reset_camera() 
        self.plotter.camera.zoom(1.3)
        self.plotter.render()

    def update_plots(self):
        t_in = np.clip(self.current_inputs['temp'], 20, 35)
        q_in = np.clip(self.current_inputs['flow'], 0, 8)
        
        rom_start_time = time.perf_counter()
        a_T = np.array([f(t_in, q_in) for f in self.rbf_T]); a_U = np.array([f(t_in, q_in) for f in self.rbf_U])
        T_recon = self.Phi_T @ a_T; U_recon = (self.Phi_U @ a_U).reshape(3, self.N).T
        Vmag = np.linalg.norm(U_recon, axis=1)
        rom_elapsed_time = time.perf_counter() - rom_start_time
        print(f"ROM calculation time: {rom_elapsed_time:.6f} seconds")

        # --- ตรวจสอบการเปลี่ยนโหมด เพื่อล้างหน้าจอเฉพาะตอนเปลี่ยนมุมมอง ---
        mode_changed = not hasattr(self, '_last_mode') or self._last_mode != self.display_mode
        submode_changed = not hasattr(self, '_last_submode') or self._last_submode != self.vector_submode
        need_rebuild = mode_changed or submode_changed

        title_str = self.display_mode
        if title_str == "Vectors": title_str = "Velocity"
        if title_str == "Do": title_str = "DO"
            
        sbar_args = {
            "title": f"{title_str}",
            "vertical": True,        
            "position_x": 0.05,      
            "position_y": 0.15,      
            "height": 0.7,           
            "width": 0.08,           
            "color": "black",        
            "title_font_size": 20,
            "label_font_size": 16,
            "fmt": "%.3f"            
        }
        
        # ล้างหน้าจอเฉพาะตอนกดปุ่มเปลี่ยนมุมมองเท่านั้น
        if need_rebuild:
            self.plotter.clear()
            self._last_mode = self.display_mode
            self._last_submode = self.vector_submode

        if self.display_mode == "Vectors":
            self.cloud.point_data['Velocity'] = Vmag
            b = self.cloud.bounds; thick = 0.005
            if self.vector_submode == "Top":
                self.cloud.point_data['vec'] = U_recon * [1, 0, 1] 
                sl_main = self.cloud.clip_box([b[0], b[1], b[3]-thick, b[3]+thick, b[4],b[5]], invert=False)
                new_mesh = sl_main.glyph(orient='vec', scale='Velocity', factor=0.15)
            else:
                self.cloud.point_data['vec'] = U_recon * [1, 1, 0]
                sl_side = self.cloud.clip_box([b[0], b[1], b[2], b[3], -thick, thick], invert=False)
                new_mesh = sl_side.glyph(orient='vec', scale='Velocity', factor=0.15) if sl_side.n_points > 0 else None

            # --- แก้ไขจุดที่กระพริบ: ใช้ name='vector_arrows' เพื่ออัปเดตลูกศรแบบ In-place ทับไปเลย ---
            if new_mesh is not None:
                self.plotter.add_mesh(new_mesh, name='vector_arrows', scalars='Velocity', cmap='turbo', scalar_bar_args=sbar_args, reset_camera=False)

        elif self.display_mode == "Probe":
            if need_rebuild:
                for name, data in self.boundaries.items():
                    self.plotter.add_mesh(data['mesh'], style='wireframe', color='black', opacity=0.15, reset_camera=False)
            
            # --- อัปเดตลูกบอลแบบ In-place ด้วย parameter name เช่นเดียวกัน ---
            if self.pt is not None:
                self.plotter.add_mesh(pv.Sphere(radius=0.005, center=self.pt), name='probe_sphere', color='red', reset_camera=False)

        else:
            s_name = self.display_mode
            data_arr = T_recon if s_name == "Temperature" else Vmag
            if s_name == "Do": data_arr = self.solver_do.predict_array(Vmag, T_recon)
            cmap = "jet" if s_name == "Temperature" else "turbo" if s_name == "Velocity" else "viridis"
            clim=self.clim_temp if s_name == "Temperature" else self.clim_vel if s_name == "Velocity" else self.clim_Do
            
            if need_rebuild:
                # สร้างวัตถุ 3D ขึ้นมาในครั้งแรก
                for name, data in self.boundaries.items():
                    data['mesh'].point_data[s_name] = data_arr[data['indices']]
                    self.plotter.add_mesh(data['mesh'], scalars=s_name, cmap=cmap, clim=clim, show_edges=False, reset_camera=False, scalar_bar_args=sbar_args)
            else:
                # ครั้งต่อไป อัปเดตแค่ค่าสี (Array) ลงในเนื้อเดิม
                for name, data in self.boundaries.items():
                    data['mesh'].point_data[s_name] = data_arr[data['indices']]
        
        # Calculation at Probe
        _, idx = self.rom_tree.query(self.pt)
        T_p, V_p = T_recon[idx], Vmag[idx]
        DO_p = self.solver_do.predict_array(np.array([V_p]), np.array([T_p]))[0]
        self.card_res_t.setText(f"{T_p:.2f}")
        self.card_res_v.setText(f"{V_p:.3f}")
        self.card_res_do.setText(f"{DO_p:.2f}")
        self.card_do_in.setText(f"{DO_p:.2f}")
        
        # ------------------------------------------
        # ตรวจสอบและอัปเดตระบบแจ้งเตือน (Alarms)
        # ------------------------------------------
        self.update_alarms(self.current_inputs['temp'], self.current_inputs['flow'], self.current_inputs['ph'], DO_p)
        
        self.plotter.render()

    def manual_probe_sync(self):
        try:
            self.pt = np.array([float(self.in_x.text()), float(self.in_y.text()), float(self.in_z.text())])
            self.update_plots()
        except: pass

    def set_view(self, mode): 
        self.display_mode = mode
        self.vector_panel.setVisible(mode == "Vectors")
        self.update_plots()
        
    def set_vector_view(self, sub): 
        self.vector_submode = sub
        self.update_plots()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = ROMApp()
    window.showMaximized()
    sys.exit(app.exec_())