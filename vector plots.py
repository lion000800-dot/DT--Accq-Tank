import pandas as pd
import numpy as np
import pyvista as pv

# --- 1. DATA LOADING ---
def load_ansys_csv(path, cols):
    with open(path, 'r') as f:
        lines = f.readlines()
    start = next(i for i, line in enumerate(lines) if "Node Num" in line)
    df = pd.read_csv(path, skiprows=start)
    df.columns = [c.strip() for c in df.columns]
    return df[cols].apply(pd.to_numeric, errors='coerce').dropna()

# Load Data
df_geo = load_ansys_csv('exportgeometry.csv', ['Node Number', 'X [ m ]', 'Y [ m ]', 'Z [ m ]'])
df_phys = load_ansys_csv('export1.csv', 
    ['Node Number', 'Temperature [ C ]', 'Velocity [ m s^-1 ]', 
     'Velocity u [ m s^-1 ]', 'Velocity v [ m s^-1 ]', 'Velocity w [ m s^-1 ]'])

merged = pd.merge(df_geo, df_phys, on='Node Number').sort_values('Node Number')

# --- 2. POINT CLOUD & SYMMETRY ---
points = merged[['X [ m ]', 'Y [ m ]', 'Z [ m ]']].values
cloud = pv.PolyData(points)
cloud.point_data['Velocity'] = merged['Velocity [ m s^-1 ]'].values
v_stack = np.column_stack((merged['Velocity u [ m s^-1 ]'].values,
                           merged['Velocity v [ m s^-1 ]'].values,
                           merged['Velocity w [ m s^-1 ]'].values))
cloud.point_data['V_Vectors'] = v_stack

mirrored_cloud = cloud.reflect((0, 0, 1), point=(0, 0, 0))
m_vectors = mirrored_cloud.point_data['V_Vectors'].copy()
m_vectors[:, 2] *= -1  
mirrored_cloud.point_data['V_Vectors'] = m_vectors

full_cloud = cloud + mirrored_cloud

# --- 3. VECTOR PROJECTIONS ---
fv = full_cloud.point_data['V_Vectors']
full_cloud.point_data['vec_xy'] = np.column_stack((fv[:, 0], fv[:, 1], np.zeros_like(fv[:, 0])))
full_cloud.point_data['vec_xz'] = np.column_stack((fv[:, 0], np.zeros_like(fv[:, 0]), fv[:, 2]))

# --- 4. VISUALIZATION ---
plotter = pv.Plotter(shape=(1, 2), window_size=[1600, 800])
plotter.set_background("white")

def add_vector_subplot(plane_normal, vec_name, view_func, title, location='center'):
    b = full_cloud.bounds
    thickness = 0.01 
    
    if plane_normal == 'z':
        # XY Plane remains at symmetry center (Z=0)
        slice_obj = full_cloud.clip_box(bounds=[b[0], b[1], b[2], b[3], -thickness, thickness], invert=False)
    else:
        # XZ Plane location logic
        if location == 'top':
            # Use the Maximum Y bound (b[3]) and subtract a small offset so we don't hit empty space
            target_y = b[3] - thickness
            slice_obj = full_cloud.clip_box(bounds=[b[0], b[1], target_y - thickness, target_y + thickness, b[4], b[5]], invert=False)
        else:
            # Default center (Y=0)
            slice_obj = full_cloud.clip_box(bounds=[b[0], b[1], -thickness, thickness, b[4], b[5]], invert=False)
    
    if slice_obj.n_points > 0:
        surface = slice_obj.delaunay_2d() 
        boundaries = surface.extract_feature_edges(boundary_edges=True, 

                                                   non_manifold_edges=False, 

                                                   feature_edges=False, 

                                                   manifold_edges=False)
        plotter.add_mesh(boundaries, color='black', line_width=2)
        
        # Color bar settings for independent subplots
        sargs = dict(title="Velocity (m/s)", position_x=0.25, position_y=0.05, vertical=False, height=0.05, width=0.5, color='black')
        
        arrows = slice_obj.glyph(orient=vec_name, scale='Velocity', factor=0.2, tolerance=0.001)
        plotter.add_mesh(arrows, scalars='Velocity', cmap='turbo', lighting=True, scalar_bar_args=sargs)
        
        view_func()
        plotter.add_text(title, color='black', font_size=12)

# Subplot 0: XY Plane (Center)
plotter.subplot(0, 0)
add_vector_subplot('z', 'vec_xy', plotter.view_xy, "XY Plane: Symmetry Center")

# Subplot 1: XZ Plane (Top)
plotter.subplot(0, 1)
# Added 'top' argument here to shift the slice location
add_vector_subplot('y', 'vec_xz', plotter.view_xz, "XZ Plane: Top Boundary Section", location='top')

plotter.show()