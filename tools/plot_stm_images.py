#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Author  :   Emma Jurmand, Mariana Rossi 
@Contact :   mariana.rossi@gmail.com
@Desc    :   A useful script for transforming cube files in STM images
'''
# author: Emma Jurmand, 2025
# based on scripts from Mariana Rossi and Dmitrii Maksimov
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from skimage.measure import marching_cubes
from scipy.stats import gaussian_kde
import argparse

def load_data(file_name):
    """load data"""
    data = np.loadtxt(file_name, skiprows=1)
    x_data = data[:, 0]
    y_data = data[:, 1]
    z_data = data[:, 2]
    v_data = data[:, 3]
    return x_data, y_data, z_data, v_data

def process_isosurface(x_data, y_data, z_data, v_data, iso_value, z_min, z_max):
    x = np.unique(x_data)
    y = np.unique(y_data)
    z = np.unique(z_data)
    V = v_data.reshape(len(x), len(y), len(z))

    verts, faces, normals, values = marching_cubes(V, level=iso_value, spacing=(x[1] - x[0], y[1] - y[0], z[1] - z[0]))
    
    iso_x, iso_y, iso_z = verts[:, 0], verts[:, 1], verts[:, 2]
    
    filtered_indices = (iso_z >= z_min) & (iso_z <= z_max)
    iso_x, iso_y, iso_z = iso_x[filtered_indices], iso_y[filtered_indices], iso_z[filtered_indices]
    
    return iso_x, iso_y, iso_z

def plot_scatter(iso_x, iso_y, iso_z):
    """interactive scatter plot"""
    x_min, x_max = iso_x.min(), iso_x.max()
    y_min, y_max = iso_y.min(), iso_y.max()
    
    plt.figure(figsize=(16, 12))
    plt.scatter(iso_x, iso_y, c=iso_z, cmap='viridis', s=5)
    plt.colorbar(label='Z Coordinate')
    plt.title('Isosurface points colored by Z')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('scaled')
    plt.show(block=False)
    
    print(f"X range: [{x_min}, {x_max}]")
    print(f"Y range: [{y_min}, {y_max}]")
    
    return x_min, x_max, y_min, y_max

def apply_periodic_transformations(iso_x, iso_y, iso_z, x_min, x_max, y_min, y_max):
    """user-specified periodic transformations."""
    x_threshold = float(input("Enter X threshold for periodic transformation (or press Enter to skip): ") or x_min)
    y_threshold = float(input("Enter Y threshold for periodic transformation (or press Enter to skip): ") or y_min)
    z_purple_max = float(input("Enter Z value for purple coloration (or press Enter to skip): ") or 7)
    
    iso_y_transformed = np.where(iso_y < y_threshold, iso_y + (y_max - y_min), iso_y)
    iso_x_transformed = np.where(iso_x < x_threshold, iso_x + (x_max - x_min), iso_x)
    
    return iso_x_transformed, iso_y_transformed, iso_z, z_purple_max

def plot_scatter_with_color(extended_x, extended_y, extended_z, z_purple_max, aspect_ratio, output_file):
    """scatter plot with purple points for specified Z range"""
    colors = np.full_like(extended_z, np.nan, dtype=object)
    colors[(extended_z >= 0) & (extended_z <= z_purple_max)] = 'violet'
    colors[(extended_z < 0) | (extended_z > z_purple_max)] = extended_z[(extended_z < 0) | (extended_z > z_purple_max)]
    
    numerical_mask = (colors != 'violet')
    numerical_colors = extended_z[numerical_mask]
    categorical_colors = (colors == 'violet')
    
    plt.figure(figsize=(16, 16/aspect_ratio))
    plt.scatter(extended_x[categorical_colors], extended_y[categorical_colors], color='violet', s=5, label=f'Z in [0, {z_purple_max}]')
    scatter = plt.scatter(extended_x[numerical_mask], extended_y[numerical_mask], c=numerical_colors, cmap='viridis', s=5)
    
    plt.colorbar(scatter, label='Z Coordinate')
    plt.title(f'Isosurface (Violet for 0 <= Z <= {z_purple_max})')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('scaled')
    plt.legend()
    plt.savefig(output_file)
    plt.close()

def plot_density(extended_x, extended_y, extended_z, output_file):
    """density plot with custom colormap"""
    xyz = np.vstack([extended_x, extended_y])
    weights = extended_z
    kde = gaussian_kde(xyz, weights=weights, bw_method=0.2)
    
    x_grid = np.linspace(extended_x.min(), extended_x.max(), 500)
    y_grid = np.linspace(extended_y.min(), extended_y.max(), 500)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)
    grid_points = np.vstack([x_mesh.ravel(), y_mesh.ravel()])
    
    density = kde(grid_points).reshape(x_mesh.shape)
    
    colors = [(0, 0, 0.1), (0, 0.5, 0), (0.6, 0.9, 0.2), (1, 1, 0)]
    custom_cmap = plt.cm.colors.LinearSegmentedColormap.from_list("custom", colors, N=100)
    
    plt.figure(figsize=(16, 16), facecolor='white')
    ax = plt.gca()
    ax.set_facecolor((0, 0, 0.2))
    
    density_threshold = 0.0024
    density_masked = np.ma.masked_where(density < density_threshold, density)
    plt.pcolormesh(x_mesh, y_mesh, density_masked, cmap=custom_cmap, shading='gouraud')
    plt.colorbar(label='Density')
    
    plt.title('Density plot')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('scaled')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Interactive STM data analysis script')
    parser.add_argument('input_file', help='input data file')
    parser.add_argument('--iso-value', type=float, default=-3e-5, help='isosurface value')
    parser.add_argument('--z-min', type=float, default=2, help='Minimum z-coordinate')
    parser.add_argument('--z-max', type=float, default=40, help='Maximum z-coordinate')
    parser.add_argument('--scatter-output', default='scatter_violet.pdf', help='output file for scatter plot')
    parser.add_argument('--density-output', default='smooth_denity.pdf', help='output file for density plot')
    
    args = parser.parse_args()
    
    x_data, y_data, z_data, v_data = load_data(args.input_file)
    
    iso_x, iso_y, iso_z = process_isosurface(
        x_data, y_data, z_data, v_data, 
        args.iso_value, args.z_min, args.z_max
    )
    
    x_min, x_max, y_min, y_max = plot_scatter(iso_x, iso_y, iso_z)
    
    periodic_x, periodic_y, periodic_z, z_purple_max = apply_periodic_transformations(
        iso_x, iso_y, iso_z, x_min, x_max, y_min, y_max
    )
    
    x_range = periodic_x.max() - periodic_x.min()
    y_range = periodic_y.max() - periodic_y.min()
    aspect_ratio = x_range / y_range
    
    plot_scatter_with_color(periodic_x, periodic_y, periodic_z, z_purple_max, aspect_ratio, args.scatter_output)
    plot_density(periodic_x, periodic_y, periodic_z, args.density_output)

if __name__ == "__main__":
    main()
