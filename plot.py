import numpy as np
import yt

NDIM = 4  # Define the number of dimensions
TT, RR, TH, PHI = 0, 1, 2, 3  # Define the component indices

def read_binary_file(filename):
    with open(filename, 'rb') as f:
        Nr = np.fromfile(f, dtype=np.int32, count=1)[0]
        Nth = np.fromfile(f, dtype=np.int32, count=1)[0]
        Nphi = np.fromfile(f, dtype=np.int32, count=1)[0]
        r_grid = np.fromfile(f, dtype=np.float64, count=Nr)
        theta_grid = np.fromfile(f, dtype=np.float64, count=Nth)
        phi_grid = np.fromfile(f, dtype=np.float64, count=Nphi)
        density = np.fromfile(f, dtype=np.float64, count=Nr*Nth*Nphi).reshape((Nr, Nth, Nphi))
        pressure = np.fromfile(f, dtype=np.float64, count=Nr*Nth*Nphi).reshape((Nr, Nth, Nphi))
        velocity = np.fromfile(f, dtype=np.float64, count=NDIM*Nr*Nth*Nphi).reshape((NDIM, Nr, Nth, Nphi))
    return Nr, Nth, Nphi, r_grid, theta_grid, phi_grid, density, pressure, velocity

filename = 'output_200.bin'
Nr, Nth, Nphi, r_grid, theta_grid, phi_grid, density, pressure, velocity = read_binary_file(filename)

# Ensure that r_grid, theta_grid, and phi_grid are 1D arrays
r_grid = r_grid.flatten()
theta_grid = theta_grid.flatten()
phi_grid = phi_grid.flatten()

# Define the data dictionary with units
data = {
    ('gas', 'density'): (density, 'g/cm**3'),
    ('gas', 'pressure'): (pressure, 'dyne/cm**2'),
    ('gas', 'velocity_x'): (velocity[RR], 'cm/s'),
    ('gas', 'velocity_y'): (velocity[TH], 'cm/s'),
    ('gas', 'velocity_z'): (velocity[PHI], 'cm/s')
}

# Define the bounding box (bbox)
bbox = np.array([
    [r_grid.min(), r_grid.max()],
    [theta_grid.min(), theta_grid.max()],
    [phi_grid.min(), phi_grid.max()]
])

# Ensure the data shapes match the expected grid shape
grid_shape = (Nr, Nth, Nphi)
assert density.shape == grid_shape
assert pressure.shape == grid_shape
assert velocity.shape == (NDIM, Nr, Nth, Nphi)

# Load the dataset with yt
ds = yt.load_uniform_grid(data, grid_shape, length_unit="cm", bbox=bbox)

# Create a slice plot
slc = yt.SlicePlot(ds, 'z', ('gas', 'density'))
slc.set_log(('gas', 'density'), False)  # Ensure the density is not in log scale
slc.set_cmap(('gas', 'density'), 'viridis')
slc.annotate_grids()
slc.annotate_velocity()
slc.annotate_timestamp(corner='upper_left', redshift=False)
slc.save("density_slice.png")  # Save the slice plot as an image

# Create a projection plot
prj = yt.ProjectionPlot(ds, 'z', ('gas', 'density'))
prj.set_log(('gas', 'density'), False)  # Ensure the density is not in log scale
prj.set_cmap(('gas', 'density'), 'plasma')
prj.annotate_grids()
prj.annotate_timestamp(corner='upper_left', redshift=False)
prj.save("density_projection.png")  # Save the projection plot as an image

# Create a phase plot
phs = yt.PhasePlot(ds.all_data(), ('gas', 'density'), ('gas', 'pressure'), ('gas', 'cell_mass'))
phs.set_log(('gas', 'density'), False)
phs.set_log(('gas', 'pressure'), False)
phs.set_xlabel('Density (g/cm^3)')
phs.set_ylabel('Pressure (dyne/cm^2)')
phs.set_cmap(('gas', 'cell_mass'), 'inferno')
phs.save("phase_plot.png")  # Save the phase plot as an image

# Create a profile plot
prof = yt.ProfilePlot(ds.all_data(), 'radius', 'density')
prof.set_log('radius', False)
prof.set_log('density', False)
prof.set_xlabel('Radius (cm)')
prof.set_ylabel('Density (g/cm^3)')
prof.save("density_profile.png")  # Save the profile plot as an image
