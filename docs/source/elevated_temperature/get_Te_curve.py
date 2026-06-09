#!/usr/bin/env python3
"""
Main driver for generating 1D potential scans and computing elevated-temperature (Te) curves.

Usage examples:
---------------
# (1) Generate displaced geometries for scan
python get_Te_main.py --generate --modes 35 36 --hessian hessian_matrix.pkl --geom geometry.in

# (2) Collect FHI-aims energies and create CSV files
python get_Te_main.py --collect --modes 35 36

# (3) Compute Te(Tphys) curve(s) from potential CSV files
python get_Te_main.py --te --potentials scan_mode_35.csv scan_mode_36.csv --req 1.11 --mred 1836.15 --plots
"""

import os
import re
import sys
import pickle
import argparse
import numpy as np
import scipy.optimize
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.cm as cm
import sympy as sp
from ase.io import read, write
from ase.calculators.aims import Aims

# =========================================================
# --- UNIT CONVERSIONS ---
# =========================================================
hatokelvin = 315775.02   # 1 a.u. -> K
ANGSTROM_TO_BOHR = 1.8897261
EV_TO_HARTREE = 0.036749322

# =========================================================
# --- FONT & MATPLOTLIB SETUP ---
# =========================================================


# Global aesthetics for JCP one-column figures
plt.rcParams.update({
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "lines.linewidth": 2,
    "lines.markersize": 3,
    # Force mathtext to use Times New Roman
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    # Avoid fallback warnings
    "font.serif": "Times New Roman",
    "font.sans-serif": "Times New Roman",
    "font.cursive": "Times New Roman",
    "font.fantasy": "Times New Roman",
    "font.monospace": "Times New Roman",
})


def load_modes_from_hessian(hessian_file, geometry_file):
    """Load Hessian, build dynamical matrix, diagonalize, and return Cartesian modes."""
    atoms = read(geometry_file, format="aims")
    natoms = len(atoms)
    masses = atoms.get_masses()  # amu

    with open(hessian_file, "rb") as f:
        H = pickle.load(f)  # Hessian (3N x 3N, eV/Å²)

    sqrt_masses = np.sqrt(np.repeat(masses, 3))
    M_inv_sqrt = np.diag(1.0 / sqrt_masses)

    D = M_inv_sqrt @ H @ M_inv_sqrt
    eigvals, eigvecs = np.linalg.eigh(D)

    cart_modes = []
    for i in range(len(eigvals)):
        v = eigvecs[:, i]
        q = M_inv_sqrt @ v  # back to Cartesian
        q = q.reshape((natoms, 3))
        q /= np.linalg.norm(q)
        cart_modes.append(q)

    return eigvals, np.array(cart_modes)


def generate_mode_displacements(
    geometry,
    modes,
    mode_labels=None,
    displacement_range_angstrom=0.5,
    n_points=51,
    traj_flag=False,
    output_dir=".",
    bond_displacement_mode=False,
    bond_pairs=None,
    bond_threshold=1.2,
    breathe_mode=False,
    breathe_pair=None
):
    """
    Generate displaced geometries along normal modes, single bonds, or
    a 'breathing' pattern that moves all bonds of the same kind simultaneously.
    """
    atoms = read(geometry, format="aims")
    os.makedirs(output_dir, exist_ok=True)

    # Handle displacement range
    if isinstance(displacement_range_angstrom, (list, tuple)):
        disp_min, disp_max = displacement_range_angstrom
    else:
        disp_min, disp_max = -displacement_range_angstrom, displacement_range_angstrom
    displacements = np.linspace(disp_min, disp_max, n_points)

    all_frames = {}

    if mode_labels is None:
        mode_labels = [f"mode_{i+1:02d}" for i in range(len(modes))]

    # --- STANDARD NORMAL MODE DISPLACEMENT ---
    if not bond_displacement_mode and not breathe_mode:
        for mode, label in zip(modes, mode_labels):
            frames = []
            mode_dir = os.path.join(output_dir, f"Scan_{label.capitalize()}")
            os.makedirs(mode_dir, exist_ok=True)

            for i, d in enumerate(displacements):
                atoms_disp = atoms.copy()
                atoms_disp.positions += d * mode
                step_dir = os.path.join(mode_dir, f"scan_step_{i:03d}")
                os.makedirs(step_dir, exist_ok=True)
                write(os.path.join(step_dir, "geometry.in"), atoms_disp, format="aims")

                if os.path.exists("control.in.tmp"):
                    os.system(f"cp control.in.tmp {os.path.join(step_dir, 'control.in')}")

                frames.append(atoms_disp)
                print(f"[{label}] Step {i:03d} | Δ = {d:.3f} Å written")

            all_frames[label] = frames
            if traj_flag:
                traj_file = os.path.join(output_dir, f"{label.lower()}_animation.xyz")
                write(traj_file, frames + frames[-2::-1], format="xyz")
                print(f"Trajectory written: {traj_file}")

    # --- SINGLE-BOND DISPLACEMENT MODE ---
    elif bond_displacement_mode and not breathe_mode:
        if not bond_pairs:
            raise ValueError("bond_pairs must be provided when bond_displacement_mode=True")

        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        natoms = len(atoms)
        found_bonds = []

        for i in range(natoms):
            for j in range(i + 1, natoms):
                dist = np.linalg.norm(positions[i] - positions[j])
                pair = tuple(sorted((symbols[i], symbols[j])))
                if pair in [tuple(sorted(p)) for p in bond_pairs] and dist < bond_threshold:
                    found_bonds.append((i, j))

        if not found_bonds:
            print("No bonds found matching the criteria.")
            return {}

        def get_neighbors(index, positions, symbols, cutoff=1.3):
            neigh = []
            for k in range(len(positions)):
                if k != index:
                    d = np.linalg.norm(positions[index] - positions[k])
                    if d < cutoff:
                        neigh.append(k)
            return neigh

        for i, (a, b) in enumerate(found_bonds):
            label = f"bond_{symbols[a]}{a+1}-{symbols[b]}{b+1}"
            bond_dir = os.path.join(output_dir, f"Scan_{label}")
            os.makedirs(bond_dir, exist_ok=True)
            frames = []

            if atoms[a].mass < atoms[b].mass:
                light_atom, heavy_atom = a, b
            elif atoms[b].mass < atoms[a].mass:
                light_atom, heavy_atom = b, a
            else:
                neigh_a = get_neighbors(a, positions, symbols)
                neigh_b = get_neighbors(b, positions, symbols)
                mass_a = sum(atoms[k].mass for k in neigh_a)
                mass_b = sum(atoms[k].mass for k in neigh_b)
                light_atom, heavy_atom = (b, a) if mass_a >= mass_b else (a, b)

            vec = positions[light_atom] - positions[heavy_atom]
            vec /= np.linalg.norm(vec)

            for k, d in enumerate(displacements):
                atoms_disp = atoms.copy()
                atoms_disp.positions[light_atom] += d * vec
                step_dir = os.path.join(bond_dir, f"scan_step_{k:03d}")
                os.makedirs(step_dir, exist_ok=True)
                write(os.path.join(step_dir, "geometry.in"), atoms_disp, format="aims")

                if os.path.exists("control.in.tmp"):
                    os.system(f"cp control.in.tmp {os.path.join(step_dir, 'control.in')}")
                frames.append(atoms_disp)
                print(f"[{label}] Step {k:03d} | Δ = {d:.3f} Å written")

            all_frames[label] = frames
            if traj_flag:
                traj_file = os.path.join(output_dir, f"{label}_animation.xyz")
                write(traj_file, frames + frames[-2::-1], format="xyz")
                print(f"Trajectory written: {traj_file}")

    # --- BOND-BREATHING MODE ---
    elif breathe_mode:
        if breathe_pair is None or len(breathe_pair) != 2:
            raise ValueError("--breathe requires a two-letter bond type (e.g., OH, NH).")

        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        natoms = len(atoms)
        pair_sorted = tuple(sorted(breathe_pair))
        found_bonds = []

        # --- find all matching bonds of the given type ---
        for i in range(natoms):
            for j in range(i + 1, natoms):
                dist = np.linalg.norm(positions[i] - positions[j])
                pair = tuple(sorted((symbols[i], symbols[j])))
                if pair == pair_sorted and dist < bond_threshold:
                    found_bonds.append((i, j))

        if not found_bonds:
            print(f"No {breathe_pair[0]}-{breathe_pair[1]} bonds found below {bond_threshold} Å.")
            return {}

        label = f"breathe_{breathe_pair[0]}{breathe_pair[1]}"
        mode_dir = os.path.join(output_dir, f"Scan_{label}")
        os.makedirs(mode_dir, exist_ok=True)
        frames = []

        # --- Helper: choose which atom moves (light atom) ---
        def select_light_atom(a, b):
            if atoms[a].mass < atoms[b].mass:
                return a, b  # light, heavy
            else:
                return b, a

        for k, d in enumerate(displacements):
            atoms_disp = atoms.copy()

            for (a, b) in found_bonds:
                light_atom, heavy_atom = select_light_atom(a, b)
                vec = positions[light_atom] - positions[heavy_atom]
                vec /= np.linalg.norm(vec)
                atoms_disp.positions[light_atom] += d * vec

            step_dir = os.path.join(mode_dir, f"scan_step_{k:03d}")
            os.makedirs(step_dir, exist_ok=True)
            write(os.path.join(step_dir, "geometry.in"), atoms_disp, format="aims")

            if os.path.exists("control.in.tmp"):
                os.system(f"cp control.in.tmp {os.path.join(step_dir, 'control.in')}")

            frames.append(atoms_disp)
            print(f"[{label}] Step {k:03d} | Δ = {d:.3f} Å written ({len(found_bonds)} bonds displaced)")

        all_frames[label] = frames

        if traj_flag:
            traj_file = os.path.join(output_dir, f"{label}_animation.xyz")
            write(traj_file, frames + frames[-2::-1], format="xyz")
            print(f"Trajectory written: {traj_file}")


    print("Geometry generation complete.")
    return all_frames


def extract_energy_from_aims_out(path):
    """
    Extracts the total energy (in eV) from an FHI-aims output file.

    Parameters
    ----------
    path : str
        Path to the FHI-aims output file.

    Returns
    -------
    energy : float or None
        Total energy in eV, or None if not found or file missing.
    """
    if not os.path.exists(path):
        return None

    # Regex pattern for total energy lines (matches both corrected/uncorrected)
    pattern = re.compile(
        r"Total\s+energy\s+of\s+the\s+DFT\s*/\s*Hartree-Fock\s*s\.c\.f\.\s*calculation\s*:?\s*([-]?\d+\.\d+)")
    try:
        with open(path, "r") as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    return float(match.group(1))
    except Exception:
        return None

    return None


def collect_and_write_energies(scan_dir, displacement_range, n_points, prefix="caf_sp"):
    """
    Collect energies from a single FHI-aims scan directory and write a CSV file.

    Parameters
    ----------
    scan_dir : str
        Directory that directly contains scan_step_XXX subdirectories.
    displacement_range : tuple of float
        (min_disp, max_disp) in Å used to reconstruct displacements.
    n_points : int
        Number of displacement points in the scan.
    prefix : str, optional
        Prefix of FHI-aims output files (e.g., 'caf_sp', 'map_sp', 'esp').
    """

    disp_min, disp_max = displacement_range
    displacements = np.linspace(disp_min, disp_max, n_points)
    energies, missing_steps = [], []

    print(f"Collecting from: {scan_dir}")

    for i in range(n_points):
        step_dir = os.path.join(scan_dir, f"scan_step_{i:03d}")
        out_file = os.path.join(step_dir, f"{prefix}_{i:03d}.out")

        E = extract_energy_from_aims_out(out_file)
        if E is None:
            missing_steps.append(i)
        energies.append(E)

    # Convert to arrays
    energies = np.array([e if e is not None else np.nan for e in energies])

    # Zero the energy (subtract minimum ignoring NaNs)
    if np.isfinite(energies).any():
        energies -= np.nanmin(energies)

    # Save to CSV
    csv_name = os.path.basename(os.path.normpath(scan_dir)) + ".csv"
    csv_path = os.path.join(os.path.dirname(scan_dir), csv_name)
    np.savetxt(
        csv_path,
        np.column_stack([displacements, energies]),
        header="Displacement (Å), Energy (eV)",
        delimiter=",",
        fmt="%.6f"
    )

    print(f"{csv_name} written in {os.path.dirname(scan_dir)}")

    if missing_steps:
        print(f"Missing steps: {missing_steps}")
    else:
        print("All points collected successfully.")



def plot_Te_curves(X_LIST, Y_LIST, LABELS=None, COLORS=None, 
                   filename='Te_curves.pdf', figsize=(3.5, 3), 
                   xlabel=r"$T_{\mathrm{phys}}$ [K]", ylabel=r"$T_e$ [K]",
                   xlim=(0, 650), ylim=(100, 650)):
    """
    Plot multiple T_e(T_phys) curves on a single figure.

    Parameters
    ----------
    X_LIST : list of np.ndarray
        List of temperature arrays (T_phys) in K.
    Y_LIST : list of np.ndarray
        List of corresponding elevated temperature arrays (T_e) in K.
    LABELS : list of str, optional
        Labels for each curve (for the legend).
    COLORS : list of color specs, optional
        Custom colors. If None, colors are taken from 'viridis' colormap.
    filename : str, optional
        Output filename for saving the figure.
    figsize : tuple, optional
        Figure size in inches (default = (3.5, 3)).
    xlabel, ylabel : str, optional
        Axis labels.
    xlim, ylim : tuple, optional
        Axis limits for x and y.
    """

    # --- Figure setup ---
    fig, ax = plt.subplots(figsize=figsize)
    cmap = plt.get_cmap('viridis')
    n_curves = len(X_LIST)

    if COLORS is None:
        start, end = 0.15, 0.85  # avoid extremes for better visual contrast
        COLORS = [cmap(start + (end - start) * i / max(1, n_curves - 1)) for i in range(n_curves)]
    COLORS = [cmap(0.2),cmap(0.6)] if n_curves==2 else COLORS  # specific for 3 curves

    if LABELS is None:
        LABELS = [f"Curve {i+1}" for i in range(n_curves)]

    # --- Plot each curve ---
    for X, Y, color, label in zip(X_LIST, Y_LIST, COLORS, LABELS):
        ax.plot(X, Y, '-', color=color, linewidth=2.5, label=label)

    # --- Axes and style ---
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.minorticks_on()
    ax.grid(which='minor', linestyle=':', alpha=0.3)
    ax.legend(frameon=False, fontsize=11)
    plt.tight_layout()
    plt.savefig(filename, format='pdf', bbox_inches='tight', dpi=600)
    plt.show()


def smooth_potential(potential_file, req, plots=False):
    # --- Load DFT scan data ---
    data = np.loadtxt(potential_file, delimiter=',', skiprows=1)  # ignore header
    disp_A, energy_eV = data.T

    # Convert to atomic units
    disp_bohr = disp_A * ANGSTROM_TO_BOHR + req
    energy_H  = energy_eV * EV_TO_HARTREE
    energy_H  = energy_H - energy_H.min()   # shift minimum to zero

    # --- Interpolated potential and derivatives ---
    potentialinterp = scipy.interpolate.CubicSpline(disp_bohr, energy_H,extrapolate=True)
    dVdxnumerical   = potentialinterp.derivative(1)
    d2Vdx2numerical = potentialinterp.derivative(2)

    # Plot the potential and interpolation
    if plots:
        xfine = np.linspace(disp_bohr.min(), disp_bohr.max(), 400)
        plt.figure(figsize=(6,4))
        plt.plot(disp_bohr, energy_H, "o", label="DFT data")
        plt.plot(xfine, potentialinterp(xfine), "-", label="Cubic spline")
        plt.xlabel("Bond length [$a_0$]")
        plt.ylabel("Energy [Hartree]")
        plt.legend()
        plt.title(f"Potential: {potential_file}")
        plt.grid(True, color="#d3d3d3")
        plt.savefig(potential_file.replace('.csv','_potential.png'), dpi=600)
        plt.show()

    return(disp_bohr,potentialinterp,dVdxnumerical,d2Vdx2numerical)

def get_rc_of_T(ktgrid,disp_bohr,dVdxnumerical,d2Vdx2numerical,mred,plots=False,label=None):
    # ---- locate the minimum:  solve dV/dx = 0 inside the data range ----
    x_min = scipy.optimize.brentq(dVdxnumerical, disp_bohr[0], disp_bohr[-1])   # refine with a root-finder

    # ---- curvature (second derivative) at the minimum ----
    kappa = d2Vdx2numerical(x_min)                  # V″(x_min)  in Hartree · bohr⁻²
    wharmnumerical = np.sqrt(kappa/mred)

    print("Vibrational frequency: ", wharmnumerical*219474.63, " cm^-1")

    w1_grid = 2 * ktgrid * np.pi              # frequency-like variable

    # Start root tracking from the minimum
    rc_guess = x_min
    solutions = []

    for w1 in w1_grid:
        def equation(rc):
            return rc + (1 / (mred * w1**2)) * dVdxnumerical(rc)
        try:
            # search in a window around the previous solution
            rc_solution = scipy.optimize.brentq(equation, rc_guess-0.5, rc_guess+0.5,
                                                xtol=1e-12, rtol=1e-10, maxiter=200)
        except ValueError:
            # fallback to fsolve near the last root
            rc_solution = scipy.optimize.fsolve(equation, rc_guess)[0]
        solutions.append(rc_solution)
        rc_guess = rc_solution   # update for continuity

    rc_array = np.array(solutions)


    # --- Extrapolated onset value (low-T plateau) ---
    rc_onset_numerical = np.mean(rc_array[:5])   # average of first 5 values for stability
    print(f"Extrapolated rc (low-T onset) = {rc_onset_numerical:.6f} a0")


    rc_interp_numerical = scipy.interpolate.interp1d(ktgrid, rc_array, kind='cubic', fill_value="extrapolate")
    kt_fine = np.logspace(np.log10(ktgrid.min()), np.log10(ktgrid.max()), 300)
    rc_fine_numerical = rc_interp_numerical(kt_fine)

    # --- Plot rc(T) if requested ---
    if plots:
        plt.figure(figsize=(6,4))
        plt.xscale('log')
        plt.plot(ktgrid*hatokelvin, rc_array, "-o", markersize=2)
        plt.plot(kt_fine*hatokelvin, rc_fine_numerical, '-', label='Interpolation',color='gold')
        plt.xlabel(r"$T$ [K]")
        plt.ylabel(r"$r_c$ [$a_0$]")
        plt.title(f"$r_c(T)$ curve: {label}")
        plt.legend()
        plt.grid(True, color="#d3d3d3")
        safe_label = label.replace('.csv', '').replace('.', '_')
        plt.savefig(f"{safe_label}_rc_of_T.png", dpi=600)
        plt.show()

    return rc_fine_numerical,wharmnumerical,rc_onset_numerical

def get_T_of_rc(ktgrid,rc_fine_numerical,plots=False,label=None):
    kt_fine = np.logspace(np.log10(ktgrid.min()), np.log10(ktgrid.max()), 300)
    
    kt_of_rc_numerical = scipy.interpolate.interp1d(
    rc_fine_numerical, kt_fine, kind='cubic', fill_value='extrapolate')

    ktnew = kt_of_rc_numerical(rc_fine_numerical)
    if plots:
        # --- Plot T as a function of rc ---
        plt.figure(figsize=(6,4))
        plt.plot(rc_fine_numerical, ktnew*hatokelvin, 'o', color='teal', label='Interpolation')
        plt.xlabel(r"$r_c$ [$a_0$]")
        plt.ylabel(r"$T$ [K]")
        plt.yscale('log')
        plt.legend()
        plt.grid(True, color="#d3d3d3")
        plt.title(f"$T(r_c)$ curve: {label}")
        safe_label = label.replace('.csv', '').replace('.', '_')
        plt.savefig(f"{safe_label}_T_of_rc.png", dpi=600)
        plt.show()
    return(kt_of_rc_numerical)

def get_sigma_of_T(dVdxnumerical,req,mred,wharm,ktgrid,plots=False):
    req = req*ANGSTROM_TO_BOHR
    # Normalized distribution
    def standatd_dev(kt, wharm, mred):
        sigma = np.sqrt(kt/(mred * wharm**2))
        return sigma

    wharm = np.sqrt(dVdxnumerical(req)/mred)
    sigma_array = 1.0 / np.sqrt(ktgrid * mred * wharm**2)

    plt.figure(figsize=(6,4))
    plt.plot(ktgrid, sigma_array, '-', color='gold')
    plt.xscale('log')
    plt.ylabel(r"$T$ (K)")
    plt.ylabel('$\\sigma$ ($a_0$)')
    plt.title('Centroid distribution width')
    plt.grid(True, color='#d3d3d3')
    plt.savefig('sigma_of_T.png', dpi=600)
    plt.show()

def get_Te_curve(potential_file, req, mred, plots=False, tmin=3.0, tmax=1000.0):
    """Compute Te(Tphys) curve and return (Tphys, Te), harmonic frequency, Tmin."""
    req = req * ANGSTROM_TO_BOHR
    disp_bohr, potentialinterp, dVdxnumerical, d2Vdx2numerical = smooth_potential(
        potential_file, req, plots=plots
    )

    # --- Temperature grid (converted from Kelvin to a.u.) ---
    kt_min = tmin / hatokelvin   # k_B T_min in a.u.
    kt_max = tmax / hatokelvin   # k_B T_max in a.u.
    ktgrid = np.logspace(np.log10(kt_min), np.log10(kt_max), num=200)

    # --- Get rc(T), frequency, and onset ---
    rc_fine_numerical, wharmnumerical, rc_onset = get_rc_of_T(
        ktgrid, disp_bohr, dVdxnumerical, d2Vdx2numerical, mred, plots, label=potential_file
    )
    kt_of_rc_numerical = get_T_of_rc(ktgrid, rc_fine_numerical, plots, label=potential_file)

    # --- Find Tmin (crossing point) ---
    def sigma(kt): return np.sqrt(kt / (mred * wharmnumerical**2))
    def delta4sigma(kt):
        r4sigma = req - 4 * sigma(kt)
        kt4sigma = kt_of_rc_numerical(r4sigma)
        return kt - kt4sigma

    deltas = np.array([delta4sigma(ktphys) for ktphys in ktgrid])
    flip = np.where(np.diff(np.signbit(deltas)))[0]
    if flip.size == 0:
        raise RuntimeError("Δ(kt) never changes sign on the supplied grid!")

    i = flip[0]
    kt_lo, kt_hi = ktgrid[i], ktgrid[i + 1]
    kt_minimum = scipy.optimize.brentq(delta4sigma, kt_lo, kt_hi,
                                       xtol=1e-12, rtol=1e-10, maxiter=100)
    Tmin_K = kt_minimum * hatokelvin

    print(f"Crossing at Tmin = {Tmin_K:.2f} K")
    print(f"Vibrational frequency = {wharmnumerical * 219474.63:.2f} cm⁻¹")

    # --- Compute Te(Tphys) curve ---
    ktphys_ktelevated = []
    for ktphys in ktgrid:
        r55sigma = req - 6 * sigma(ktphys)
        #rmin, rmax = rc_fine_numerical.min(), rc_fine_numerical.max()
        #if r55sigma < rmin or r55sigma > rmax:
        #    kt55sigma = ktphys
        #else:
        kt55sigma = kt_of_rc_numerical(r55sigma)

        if kt55sigma < kt_minimum:
            ktphys_ktelevated.append([ktphys * hatokelvin, Tmin_K])
        elif kt55sigma >= kt_minimum and kt55sigma >= ktphys:
            ktphys_ktelevated.append([ktphys * hatokelvin, kt55sigma * hatokelvin])
        elif kt55sigma < ktphys:
            ktphys_ktelevated.append([ktphys * hatokelvin, ktphys * hatokelvin])

    ktphys_ktelevated = np.array(ktphys_ktelevated)
    xtt, ytt = ktphys_ktelevated.T

    # Return Tphys, Te, frequency (cm⁻¹), Tmin (K)
    return xtt, ytt, wharmnumerical * 219474.63, Tmin_K

def generate_morse_potential(
    De=0.18748, a=1.1605, Re=1.8324, rmin=-0.5, rmax=0.5, npoints=51, plots=False
):
    """
    Generate a Morse potential and its derivatives.
    The displacement range is given in Å relative to Re (so total range = Re + Δr).

    Parameters
    ----------
    De : float
        Well depth in eV (converted internally to Hartree).
    a : float
        Range parameter in Å⁻¹ (converted internally to Bohr⁻¹).
    Re : float
        Equilibrium bond length in Å (converted internally to Bohr).
    rmin, rmax : float
        Displacement range relative to Re (Å).
    npoints : int
        Number of grid points.
    plots : bool
        If True, show potential plot.

    Returns
    -------
    disp_bohr : np.ndarray
        Grid points in Bohr.
    potentialinterp : callable
        V(r) in Hartree.
    dVdxnumerical : callable
        First derivative.
    d2Vdx2numerical : callable
        Second derivative.
    """
    # Convert to atomic units
    De_H = De * EV_TO_HARTREE
    a_B = a / ANGSTROM_TO_BOHR
    Re_B = Re * ANGSTROM_TO_BOHR
    dr_B = np.linspace(rmin, rmax, npoints) * ANGSTROM_TO_BOHR
    r_grid = Re_B + dr_B

    # Define symbolic Morse potential
    r = sp.Symbol('r')
    V = De_H * (1 - sp.exp(-a_B * (r - Re_B)))**2
    dVdr = sp.diff(V, r)
    d2Vdr2 = sp.diff(dVdr, r)

    # Lambdify
    V_num = sp.lambdify(r, V, modules='numpy')
    dVdr_num = sp.lambdify(r, dVdr, modules='numpy')
    d2Vdr2_num = sp.lambdify(r, d2Vdr2, modules='numpy')

    # Evaluate
    energy_H = V_num(r_grid)
    energy_H -= np.min(energy_H)

    if plots:
        plt.figure(figsize=(6, 4))
        plt.plot(r_grid, energy_H, '-', color='purple', label='Morse potential')
        plt.xlabel(r"$r$ [$a_0$]")
        plt.ylabel("Energy [Hartree]")
        plt.title("Analytical Morse potential")
        plt.legend()
        plt.grid(True, color="#d3d3d3")
        plt.savefig("morse_potential.png", dpi=600)
        plt.show()

    potentialinterp = lambda x: V_num(x)
    return r_grid, potentialinterp, dVdr_num, d2Vdr2_num

def get_Te_curve_from_morse(
    disp_bohr, dVdxnumerical, d2Vdx2numerical, req, mred,
    plots=False, tmin=3.0, tmax=1000.0, label="Morse"
):
    """Compute Te(Tphys) for an analytical Morse potential already in memory."""
    req = req * ANGSTROM_TO_BOHR

    # Temperature grid (converted from Kelvin to a.u.)
    kt_min = tmin / hatokelvin
    kt_max = tmax / hatokelvin
    ktgrid = np.logspace(np.log10(kt_min), np.log10(kt_max), num=200)

    rc_fine_numerical, wharmnumerical, rc_onset = get_rc_of_T(
        ktgrid, disp_bohr, dVdxnumerical, d2Vdx2numerical, mred, plots, label=label
    )
    kt_of_rc_numerical = get_T_of_rc(ktgrid, rc_fine_numerical, plots, label=label)

    # Find Tmin
    def sigma(kt): return np.sqrt(kt / (mred * wharmnumerical**2))
    def delta4sigma(kt):
        r4sigma = req - 4 * sigma(kt)
        kt4sigma = kt_of_rc_numerical(r4sigma)
        return kt - kt4sigma

    deltas = np.array([delta4sigma(ktphys) for ktphys in ktgrid])
    flip = np.where(np.diff(np.signbit(deltas)))[0]
    if flip.size == 0:
        raise RuntimeError("Δ(kt) never changes sign for Morse potential!")

    i = flip[0]
    kt_minimum = scipy.optimize.brentq(delta4sigma, ktgrid[i], ktgrid[i + 1])
    Tmin_K = kt_minimum * hatokelvin

    print(f"Crossing at Tmin = {Tmin_K:.2f} K")
    print(f"Vibrational frequency = {wharmnumerical * 219474.63:.2f} cm⁻¹")

    # Compute Te(Tphys)
    ktphys_ktelevated = []
    for ktphys in ktgrid:
        r55sigma = req - 6 * sigma(ktphys)
        kt55sigma = kt_of_rc_numerical(r55sigma)
        if kt55sigma < kt_minimum:
            ktphys_ktelevated.append([ktphys * hatokelvin, Tmin_K])
        elif kt55sigma >= kt_minimum and kt55sigma >= ktphys:
            ktphys_ktelevated.append([ktphys * hatokelvin, kt55sigma * hatokelvin])
        elif kt55sigma < ktphys:
            ktphys_ktelevated.append([ktphys * hatokelvin, ktphys * hatokelvin])

    xtt, ytt = np.array(ktphys_ktelevated).T
    return xtt, ytt, wharmnumerical * 219474.63, Tmin_K


# =========================================================
# --- MAIN ROUTINE ---
# =========================================================
def main():
    parser = argparse.ArgumentParser(
        description="Perform 1D mode scans and compute Te(Tphys) curves."
    )
    parser.add_argument("--generate", action="store_true",
                        help="Generate displaced geometries along normal modes.")
    parser.add_argument("--collect", action="store_true",
                        help="Collect FHI-aims energies and write CSV potential files.")
    parser.add_argument("--te", action="store_true",
                        help="Compute Te(Tphys) curve(s) from potential scan CSV files.")
    parser.add_argument("--modes", nargs="+", type=int, default=None,
                        help="Indices of normal modes to process (e.g., 35 36).")
    parser.add_argument("--hessian", type=str, 
                        help="Path to Hessian pickle file.")
    parser.add_argument("--geom", type=str, 
                        help="Path to geometry.in file.")
    parser.add_argument("--potentials", nargs="+", default=None,
                        help="List of CSV potential files for --te mode.")
    parser.add_argument("--req", type=float, default=None,
                        help="Equilibrium bond length (Å).")
    parser.add_argument("--mred", type=float, default=None,
                        help="Reduced mass in atomic units.")
    parser.add_argument("--req-list", nargs="+", type=float,
                        help="List of equilibrium bond lengths (Å) for multiple potentials.")
    parser.add_argument("--mred-list", nargs="+", type=float,
                        help="List of reduced masses (a.u.) for multiple potentials.")
    parser.add_argument("--labels", nargs="+", type=str,
                    help="Custom labels for each potential, including Morse if present. "
                         "Example: --labels Morse CAF_Bond1 CAF_Bond2")
    parser.add_argument("--plots", action="store_true",
                        help="Enable diagnostic plots.")
    parser.add_argument("--output", type=str, default=".",
                        help="Output directory for generated data.")
    parser.add_argument("--traj", action="store_true",
                        help="Generate animation trajectories (.xyz) for each mode during geometry generation.")
    parser.add_argument("--prefix", type=str, default="aims",
                        help="Prefix for FHI-aims output files (default: aims).")
    parser.add_argument("--input", type=str, default=".",
                        help="Input directory containing folders with FHI-aims outputs.")
    
    # Displacement grid settings
    parser.add_argument("--range", nargs=2, type=float, default=[-0.5, 0.5],
                        help="Displacement range in Å (min max). Default: -0.5 0.5")
    parser.add_argument("--npoints", type=int, default=51,
                        help="Number of displacement points. Default: 51")
    
    # Bond-based generation
    parser.add_argument("--bonds", nargs="+", type=str, 
                        help="List of bonds to scan (e.g., OH NH CH)")
    parser.add_argument("--bond-threshold", type=float, default=1.2, 
                        help="Distance cutoff for bond detection (Å)")
    
    #Temperature range for Te curves
    parser.add_argument("--tmin", type=float, default=3.0, 
                        help="Minimum physical temperature (K) for Te analysis. Default: 3 K")
    parser.add_argument("--tmax", type=float, default=1000.0,
                        help="Maximum physical temperature (K) for Te analysis. Default: 1000 K")
    parser.add_argument("--tphys", type=float, default=300,
                        help="Physical temperature of the system. Default: 300 K")

    parser.add_argument("--breathe", type=str,
                    help="Apply symmetric displacement to all bonds of the given type (e.g., OH, NH).")

    parser.add_argument(
        "--morse",
        nargs="*",
        type=float,
        metavar="params",
        help=("Use an analytical Morse potential (atomic units). "
            "Defaults: De=0.18748 Ha, a=1.1605 1/Bohr, Re=1.8324 Bohr, "
            "rmin=-0.5 Bohr, rmax=0.5 Bohr, npoints=51. "
            "Example: --morse 0.18748 1.1605 1.8324 -0.5 0.5 51")
    )


    args = parser.parse_args()

    # ========== MODE 1: GENERATE GEOMETRIES ==========
    # ===============================================================
    # === MODE-BASED DISPLACEMENT GENERATION ========================
    # ===============================================================
    if args.generate and args.modes:
        if not args.hessian:
            raise ValueError("You must specify --hessian when using --modes.")
        eigvals, cart_modes = load_modes_from_hessian(args.hessian, args.geom)
        selected_modes = [cart_modes[i - 1] for i in args.modes]  # 1-indexed
        mode_labels = [f"mode_{i}" for i in args.modes]
        generate_mode_displacements(
            geometry=args.geom,
            modes=selected_modes,
            mode_labels=mode_labels,
            displacement_range_angstrom=(args.range[0], args.range[1]),
            n_points=args.npoints,
            traj_flag=args.traj,
            output_dir=args.output,
        )

    # ===============================================================
    # === BOND-BASED DISPLACEMENT GENERATION ========================
    # ===============================================================
    elif args.generate and args.breathe:
        pair = args.breathe.strip().upper()
        if len(pair) != 2:
            raise ValueError("--breathe must be a two-letter bond identifier like OH or NH.")
        generate_mode_displacements(
            geometry=args.geom,
            modes=[],
            displacement_range_angstrom=(args.range[0], args.range[1]),
            n_points=args.npoints,
            traj_flag=args.traj,
            output_dir=args.output,
            bond_displacement_mode=False,
            breathe_mode=True,
            breathe_pair=(pair[0], pair[1]),
            bond_threshold=args.bond_threshold,
        )

    elif args.generate and args.bonds:
        bond_pairs = [(b[0], b[1]) for b in args.bonds if len(b) == 2]
        generate_mode_displacements(
            geometry=args.geom,
            modes=[],  # ignored
            displacement_range_angstrom=(args.range[0], args.range[1]),
            n_points=args.npoints,
            traj_flag=args.traj,
            output_dir=args.output,
            bond_displacement_mode=True,
            bond_pairs=bond_pairs,
            bond_threshold=args.bond_threshold,
        )

    # ========== MODE 2: COLLECT ENERGIES ==========
    if args.collect:
        collect_and_write_energies(
            scan_dir=args.input,
            displacement_range=tuple(args.range) if args.range else (-0.5, 0.5),
            n_points=args.npoints,
            prefix=args.prefix
        )
        
        print(" Energy collection and CSV writing complete.")
        return

    
     # ================================================================
    # === MODE 3: COMPUTE Te CURVES ==================================
    # ================================================================
    elif args.te:
        os.makedirs(args.output, exist_ok=True)

        X_LIST, Y_LIST, LABELS = [], [], []
        freqs_list, Tmin_list = [], []

        morse_active = args.morse is not None

        # ================================================================
        # --- Handle analytical Morse potential (optional) ---
        # ================================================================
        if morse_active:
            print("\n=== Analytical Morse potential detected ===")

            # Default Morse parameters (same units as before)
            De_eV, a_Ainv, Re_A, rmin_A, rmax_A, npoints = 5.099, 2.192, 0.970, -0.5, 0.5, 51

            # Override defaults if user provided any
            if len(args.morse) > 0:
                defaults = [De_eV, a_Ainv, Re_A, rmin_A, rmax_A, npoints]
                values = args.morse + defaults[len(args.morse):]
                De_eV, a_Ainv, Re_A, rmin_A, rmax_A, npoints = values[:6]

            # Label first (before usage)
            if args.labels and len(args.labels) > 0:
                morse_label = args.labels[0]
            else:
                morse_label = f"Morse_De{De_eV:.3f}_a{a_Ainv:.3f}_Re{Re_A:.3f}"

            # Generate Morse potential
            disp_bohr, potentialinterp, dVdxnumerical, d2Vdx2numerical = generate_morse_potential(
                De_eV, a_Ainv, Re_A, rmin_A, rmax_A, int(npoints), plots=args.plots
            )

            # Equilibrium distance for Te curve interface (Å)
            req_A = Re_A

            # Detect reduced mass
            if args.mred is not None:
                mred_use = args.mred
            elif args.mred_list is not None and len(args.mred_list) > 0:
                try:
                    mred_use = float(args.mred_list[0])
                except (TypeError, ValueError):
                    mred_use = None
            else:
                mred_use = None

            if mred_use is None:
                print("Error: no valid reduced mass detected (use --mred or --mred-list).")
                sys.exit(1)

            # Compute Te(Tphys) for analytical Morse potential
            xtt, ytt, freq_cm1, Tmin_K = get_Te_curve_from_morse(
                disp_bohr, dVdxnumerical, d2Vdx2numerical,
                req=req_A, mred=mred_use,
                plots=args.plots, tmin=args.tmin, tmax=args.tmax,
                label=morse_label
            )

            X_LIST.append(xtt)
            Y_LIST.append(ytt)
            LABELS.append(morse_label)
            freqs_list.append(freq_cm1)
            Tmin_list.append(Tmin_K)

        # ================================================================
        # --- Handle CSV-based potentials (standard behavior) ---
        # ================================================================
        # Only check req/mred if no Morse potential was provided
        if not morse_active:
            if (args.req is None and args.req_list is None) or (args.mred is None and args.mred_list is None):
                print(" Please provide either --req/--mred or --req-list/--mred-list for Te computation.")
                sys.exit(1)

        # Detect CSV files
        if args.potentials is None:
            csv_files = sorted([
                os.path.join(args.input, f) for f in os.listdir(args.input)
                if f.endswith(".csv")
            ])
        else:
            csv_files = args.potentials

        # Allow Morse-only run
        if not csv_files and not morse_active:
            print(f"No .csv files found in {args.input}")
            sys.exit(1)

        # Process CSV potentials
        if csv_files:
            nfiles = len(csv_files)

            # Handle req/mred lists
            if args.req_list is not None:
                if len(args.req_list) != nfiles:
                    print(" Error: Number of req values (--req-list) must match number of potential files.")
                    sys.exit(1)
                reqs = args.req_list
            else:
                reqs = [args.req] * nfiles

            if args.mred_list is not None:
                if len(args.mred_list) != nfiles:
                    print(" Error: Number of mred values (--mred-list) must match number of potential files.")
                    sys.exit(1)
                mreds = args.mred_list
            else:
                mreds = [args.mred] * nfiles

            print("\n=== Te(Tphys) computation summary ===")
            for i, (pot_file, req, mred) in enumerate(zip(csv_files, reqs, mreds)):
                print(f" -> {os.path.basename(pot_file)}  req={req:.3f} Å,  mred={mred:.3f} a.u.")
                xtt, ytt, freq_cm1, Tmin_K = get_Te_curve(
                    pot_file, req, mred, plots=args.plots, tmin=args.tmin, tmax=args.tmax
                )
                X_LIST.append(xtt)
                Y_LIST.append(ytt)

                # Custom labels (if provided after Morse)
                label_index = i + (1 if morse_active else 0)
                if args.labels and len(args.labels) > label_index:
                    custom_label = args.labels[label_index]
                    LABELS.append(custom_label)
                else:
                    LABELS.append(os.path.basename(pot_file).replace('.csv', ''))

                freqs_list.append(freq_cm1)
                Tmin_list.append(Tmin_K)
            print("=====================================\n")

        # ================================================================
        # --- Plot and export results ---
        # ================================================================
        plot_path = os.path.join(args.output, "Te_curves.pdf")
        plot_Te_curves(X_LIST, Y_LIST, LABELS=LABELS,
                       filename=plot_path,
                       figsize=(3.5, 3), xlim=(0, 650), ylim=(100, 650))
        print(f"\nTe(Tphys) computation complete. Plot saved as {plot_path}")

        # --- Save each Te(Tphys) curve as a .txt file ---
        print("\nSaving Te(Tphys) data to text files:")
        for idx, (label, xtt, ytt) in enumerate(zip(LABELS, X_LIST, Y_LIST)):
            if label.startswith("Morse"):
                req_use = Re_A
                mred_use = mred_use
            else:
                offset = 1 if morse_active else 0
                req_use = reqs[idx - offset]
                mred_use = mreds[idx - offset]

            out_txt = os.path.join(args.output, f"{label}_Te_curve.txt")
            np.savetxt(out_txt,
                       np.column_stack([xtt, ytt]),
                       header=f"T_phys(K)   T_e(K)\nreq = {req_use:.4f} Å, mred = {mred_use:.2f} a.u.",
                       fmt="%12.6f")
            print(f"  -> {os.path.basename(out_txt)} written")

        # ================================================================
        # --- Summary table ---
        # ================================================================
        print("\n=== Summary of Te(Tphys) computations ===")
        print(f"{'Label':<30} {'req (Å)':>10} {'mred (a.u.)':>12} {'freq (cm⁻¹)':>14} "
              f"{'Tmin (K)':>10} {f'Te({args.tphys:.0f}K) (K)':>15}")
        print("-" * 96)

        summary_data = []
        for idx, (label, xtt, ytt, freq, Tmin) in enumerate(zip(LABELS, X_LIST, Y_LIST, freqs_list, Tmin_list)):
            if label.startswith("Morse"):
                req_use = Re_A
                mred_use = mred_use
            else:
                offset = 1 if morse_active else 0
                req_use = reqs[idx - offset]
                mred_use = mreds[idx - offset]

            te_interp = np.interp(args.tphys, xtt, ytt)
            summary_data.append((label, req_use, mred_use, freq, Tmin, te_interp))
            print(f"{label:<30} {req_use:>10.4f} {mred_use:>12.2f} {freq:>14.1f} "
                  f"{Tmin:>10.1f} {te_interp:>15.1f}")
        print("-" * 96)

        # --- LaTeX output ---
        latex_path = os.path.join(args.output, "Te_summary.tex")
        with open(latex_path, "w") as f:
            f.write("\\begin{table}[h!]\n\\centering\n")
            f.write("\\caption{Summary of $T_e(T_{\\mathrm{phys}})$ computations.}\n")
            f.write("\\begin{tabular}{lccccc}\n\\hline\n")
            f.write("Label & $r_{\\mathrm{eq}}$ (\\AA) & $\\mu$ (a.u.) & "
                    "$\\omega_{\\mathrm{harm}}$ (cm$^{-1}$) & $T_{\\min}$ (K) & "
                    f"$T_e({args.tphys:.0f}$ K) \\\\\n\\hline\n")
            for label, req_use, mred_use, freq, Tmin, Te_phys in summary_data:
                label_tex = label.replace('_', '\\_')
                f.write(f"{label_tex} & {req_use:.4f} & {mred_use:.2f} & "
                        f"{freq:.1f} & {Tmin:.1f} & {Te_phys:.1f} \\\\\n")
            f.write("\\hline\n\\end{tabular}\n\\end{table}\n")
        print(f"LaTeX table written to: {latex_path}")


# =========================================================
# --- ENTRY POINT ---
# =========================================================
if __name__ == "__main__":
    main()
