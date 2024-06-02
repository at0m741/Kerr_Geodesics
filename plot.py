import yt
import numpy as np
import h5py
import os

def spherical_to_cartesian_velocity(vr, vtheta, vphi, theta, phi):
    vx = vr * np.sin(theta) * np.cos(phi) + vtheta * np.cos(theta) * np.cos(phi) - vphi * np.sin(phi)
    vy = vr * np.sin(theta) * np.sin(phi) + vtheta * np.cos(theta) * np.sin(phi) + vphi * np.cos(phi)
    vz = vr * np.cos(theta) - vtheta * np.sin(theta)
    return vx, vy, vz

def load_hdf5_data(filename):
    with h5py.File(filename, 'r') as f:
        density = f['density'][:]
        pressure = f['pressure'][:]
        velocity = f['velocity'][:]
        r_grid = f['r_grid'][:]
        theta_grid = f['theta_grid'][:]
        phi_grid = f['phi_grid'][:]
    Nr = len(r_grid)
    Nth = len(theta_grid)
    Nphi = len(phi_grid)
    
    density = density.reshape((Nr, Nth, Nphi))
    pressure = pressure.reshape((Nr, Nth, Nphi))
    velocity = velocity.reshape((4, Nr, Nth, Nphi))
    
    return density, pressure, velocity, r_grid, theta_grid, phi_grid

def create_yt_dataset(density, pressure, velocity, r_grid, theta_grid, phi_grid):
    theta, phi = np.meshgrid(theta_grid, phi_grid, indexing='ij')
    theta = np.broadcast_to(theta, velocity[1].shape)
    phi = np.broadcast_to(phi, velocity[1].shape)
    
    vx, vy, vz = spherical_to_cartesian_velocity(velocity[1], velocity[2], velocity[3], theta, phi)

    data = {
        ('gas', 'density'): (density, 'g/cm**3'),
        ('gas', 'pressure'): (pressure, 'dyne/cm**2'),
        ('gas', 'velocity_x'): (vx, 'cm/s'),
        ('gas', 'velocity_y'): (vy, 'cm/s'),
        ('gas', 'velocity_z'): (vz, 'cm/s'),
    }
    
    bbox = np.array([[r_grid.min(), r_grid.max()],
                     [theta_grid.min(), theta_grid.max()],
                     [phi_grid.min(), phi_grid.max()]])
    
    ds = yt.load_uniform_grid(data, density.shape, length_unit="pc", bbox=bbox)
    return ds

def generate_images(input_files, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, filename in enumerate(input_files):
        density, pressure, velocity, r_grid, theta_grid, phi_grid = load_hdf5_data(filename)
        ds = create_yt_dataset(density, pressure, velocity, r_grid, theta_grid, phi_grid)

        slc = yt.SlicePlot(ds, 'z', ('gas', 'density'))
        slc.set_log('density', False)
        slc.set_cmap(('gas', 'density'), 'inferno')
        slc.set_width((r_grid.max() - r_grid.min(), 'pc'))
        width = 1086 
        aspect_ratio = (theta_grid.max() - theta_grid.min()) / (r_grid.max() - r_grid.min())
        height = int(width * aspect_ratio)
        height += height % 3

        slc.set_buff_size((width, height))
        slc.annotate_title(f'ADAF simulation: density slice at t = {i}')
        output_filename = os.path.join(output_dir, f'density_slice_{i:04d}.png')
        slc.save(output_filename)

input_files = [f'output_{i}.h5' for i in range(200)]
output_dir = 'frames'
generate_images(input_files, output_dir)
