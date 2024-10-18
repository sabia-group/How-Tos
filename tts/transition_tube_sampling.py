from typing import Iterable, Sequence
import warnings

import numpy as np
import scipy.constants as cs
from scipy.interpolate import CubicSpline
from scipy.linalg import orth


class NormalModes:
    """A representation of normal modes.

    Provided a Hessian matrix, calculate normal mode vectors and corresponding frequencies,
    optionally filtering out translational and rotational degrees of freedom or unstable modes.

    Normal modes and frequencies are made available in mass-weighted atomic units.
    """

    @classmethod
    def from_cp2k(cls, fn_in: str, n_atoms: int, masses: np.ndarray = None):
        """Read a CP2K Hessian file and return a NormalModes object."""
        hessian = load_cp2k_Hessian(fn_in, n_atoms)
        return cls(hessian, masses)

    @classmethod
    def from_ase(cls, hessian_ase: np.ndarray, masses: np.ndarray):
        """Take an ASE Hessian, convert units, mass-weigh and return a NormalModes object."""
        hessian = convert_ASE_hessian(hessian_ase, masses)
        return cls(hessian, masses)

    def __init__(self, hessian: np.ndarray, masses: np.ndarray = None):
        """Store the Hessian matrix to prepare for normal modes calculation.

        To obtain modes and frequencies, `update_normal_modes()` needs to be called by hand after construction.

        Arguments:
            hessian: 3N x 3N mass-weighted cartesian Hessian matrix in atomic units
            masses: N masses, one for each atom, in atomic mass constant units (where carbon is ~12)
        """

        # store the Hessian
        self.hessian = hessian

        # store masses
        self.masses = masses

        # normal modes haven't been calculated yet
        self._modes = None

    @property
    def modes(self):
        """Normal mode vectors in atomic units."""
        if self._modes is None:
            raise Exception('Normal modes have not been computed yet.')
        return self._modes.copy()

    @property
    def frequencies(self):
        """Normal mode frequencies in atomic units."""
        if self._modes is None:
            raise Exception('Normal modes have not been computed yet.')
        return self._frequencies.copy()

    @property
    def wavenumbers(self):
        """Normal mode wavenumbers in spectroscopic units (cm^-1)."""
        if self._modes is None:
            raise Exception('Normal modes have not been computed yet.')

        # convert a.u. to SI units (Hz)
        omega = self._frequencies * cs.physical_constants['atomic unit of energy'][0] / cs.hbar

        # convert angular frequency to linear
        nu = omega / (2 * np.pi)

        # convert `nu` to wavenumbers in cm^-1
        wavenumbers = nu / cs.c
        wavenumbers *= 1e-2
        return wavenumbers

    def update_normal_modes(
        self,
        symmetrize: bool = False,
        order_SP: int = 0,
        discard_imag_freq: bool = True,
        eps_freq: float = 5.0e-4,
        n_TR: int = 6
    ) -> None:
        """Diagonalize the Hessian to get normal modes and their frequencies.

        Optionally also filter out translational and rotational degrees of freedom.

        Arguments:
            symmetrize: activate symmetrization of the Hessian
            order_SP: expected saddle point order, number of finite imaginary frequencies at the configuration
            discard_imag_freq: discard unstable modes when `order_SP` > 0
            eps_freq: numerical threshold to check frequencies in a.u., 5e-4 ~ 100 cm^-1.
            n_TR: how many modes to filter to remove rotations and translations
        """

        # use a symmetrized Hessian if requested
        if symmetrize:
            hessian = 0.5 * (self.hessian + self.hessian.T)
        else:
            hessian = self.hessian

        # diagonalize
        ev, modes = np.linalg.eigh(hessian)

        # get frequencies from `ev` (omega^2), possibly imaginary
        frequencies = np.sqrt(ev.astype('complex'))

        # check that the number of finite imaginary frequencies checks out
        n_imag = (np.imag(frequencies) > eps_freq).sum()
        if n_imag != order_SP:
            warnings.warn(f'Number of imaginary frequencies ({n_imag}) differs from the expected SP order ({order_SP})')

        # transpose modes so that they're ordered along the first axis
        modes = modes.T

        # filter modes as requested
        mask = np.ones(frequencies.size, dtype='bool')
        if discard_imag_freq:
            mask[:order_SP] = False
        mask[order_SP:order_SP+n_TR] = False
        frequencies = frequencies[mask]
        modes = modes[mask]

        # lift mass weighing from modes if masses are assigned
        if self.masses is not None:
            amu = cs.physical_constants['atomic mass constant'][0]
            mass_atomic = cs.physical_constants['atomic unit of mass'][0]
            amu_to_me = amu / mass_atomic
            masses = self.masses * amu_to_me
            modes = modes / np.sqrt(np.repeat(masses, 3))

        # recast remaining real frequencies after filtering as real numbers
        frequencies = np.real(frequencies)

        # store the settings and results
        self._symmetrize = symmetrize
        self._n_TR = n_TR
        self._order_SP = order_SP
        self._discard_imag_freq = discard_imag_freq
        self._eps_freq = eps_freq
        self._modes = modes.copy()
        self._frequencies = frequencies.copy()


def load_cp2k_Hessian(fn_in: str, n_atoms: int, correct_CP2K_factor: bool = True) -> np.ndarray:
    """Load the raw Hessian matrix from CP2K output.

    Arguments:
       fn_in: Path to the CP2K log file
       n_atoms: Number of atoms in the studied molecule
       symmetrize: Whether to symmetrize the Hessian after reading
       correct_CP2K_factor: remove CP2K redundant multiplicative 1e6 factor from the Hessian, so that it is in a.u.

    Returns:
       Hessian matrix as a NumPy array in atomic units
    """

    # variables such as stride of data, regular expressions and collectors
    ndof = 3 * n_atoms
    stride = ndof + 2
    begin_regex = 'Hessian in cartesian coordinates'
    end_regex = 'Cartesian Low frequencies'

    # read all the lines of the log file, put them in a list
    with open(fn_in) as f_in:
        lines = f_in.readlines()

    # locate the beginning of the Hessian section
    i_start = None
    for i, line in enumerate(lines):
        if begin_regex in line:
            i_start = i + 3
            break
    if i_start is None:
        raise ValueError('Start of Hessian section not found.')

    # locate the end of the Hessian section
    i_end = None
    for i, line in enumerate(lines):
        if end_regex in line:
            i_end = i - 1
            break
    if i_end is None:
        raise ValueError('End of Hessian section not found.')

    # check that the Hessian section is not empty
    if i_end <= i_start:
        raise ValueError('Hessian section is empty.')

    # read in Hessian data
    hessian = np.empty((ndof, 0))
    for i in range(i_start, i_end, stride):
        data = []
        for j in range(i, i+stride-2):
            data.extend(lines[j].split()[2:])

        data = np.array(data, dtype=float).reshape(ndof, int(len(data) / ndof))

        hessian = np.append(hessian, data, axis=1)

    # check resulting shape is as expected
    assert hessian.shape == (ndof, ndof), 'Hessian matrix does not have the right shape of (3*n_atoms, 3*n_atoms).'

    if correct_CP2K_factor:
        hessian *= 1e-6

    return hessian


def convert_ASE_hessian(hessian_ase: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Convert the ASE-generated Hessian matrix into mass-weighted atomic units.

    Arguments:
        ase_hessian: 2D symmetric hessian matrix in eV * Angstrom**-2
        masses: masses for each atom ordered as in `pos` in atomic mass constant units

    Returns:
        hessian: 2D symmetric Hessian in E_h * a_0**-2 * m_e**-1
    """

    # constants
    e_h = cs.physical_constants['atomic unit of energy'][0]
    eV = cs.physical_constants['electron volt'][0]
    a_0 = cs.physical_constants['atomic unit of length'][0]
    angstrom = 1e-10
    amu = cs.physical_constants['atomic mass constant'][0]
    m_e = cs.physical_constants['atomic unit of mass'][0]

    # convert unit
    hessian = hessian_ase
    hessian *= (eV / e_h)
    hessian /= (angstrom / a_0)**2

    # define auxiliary mass matrix m_ij = sqrt(m_i)*sqrt(m_j)
    aux_mass_matrix = np.outer(np.repeat(masses, 3), np.repeat(masses, 3))

    # mass-weigh the hessian
    hessian /= (np.sqrt(aux_mass_matrix) * amu / m_e)

    return hessian


def rho_harmonic(q: np.ndarray, frequency: float, beta: float) -> np.ndarray:
    """Calculate the harmonic thermal density.

    For a given harmonic mode with frequency `frequency` and inverse temperature `beta`,
    calculate the thermal probability density along the values of the coordinate `q`.
    This can be used in the classical case or in the quantum case with an effective inverse temperature.
    Be careful about consistent units.

    Arguments:
        q: 1D array of coordinate values
        frequency: frequency of the mode
        beta: inverse temperature

    Returns:
        1D array of the density values
    """

    return np.sqrt(0.5 * beta * frequency**2 / np.pi) * np.exp(-0.5 * beta * frequency**2 * q**2)


class NormalModeSampling:
    """Sample geometries based on normal modes.

    Generate geometry samples around a single minimum or around a minimum energy path.
    The weight of the sampling is uniform along the MEP between the endpoints.
    Perpendicular to the MEP and beyong the endpoints, the weighting is thermal, based on either
    classical or quantum thermal probability densities.
    """

    def __init__(self, positions: np.ndarray, normal_modes: dict):
        """
        Arguments:
            positions: Atomic positions of all images in atomic units (Bohr radii)
            normal_modes: Dictionary of `NormalModes` objects with keys being
                the relevant indices of images that the normal modes correspond to
        """

        # store the positions and normal modes
        self._positions = positions
        self._normal_modes = normal_modes

        # set MEP length
        self._number_of_replicas = len(self._positions)

        # check that we have at least one image and that we have enough modes
        assert self._number_of_replicas >= 1
        if self._number_of_replicas == 1:
            assert len(self._normal_modes) == 1
        else:
            assert len(self._normal_modes) >= 2

        # indices of images that have normal modes available, sorted by increasing index
        idx_nm = np.array(sorted(list(self._normal_modes.keys())), dtype=int)
        self._idx_nm = idx_nm

        # If we are sampling around an MEP, rather than just a single minimum,
        # prepapre MEP parametrization and interpolation.
        if self._number_of_replicas > 1:
            # get the physical length of mep in a0 and parametrize the MEP by dimensionless xi
            d_vec = np.diff(positions, axis=0)
            d_norm = np.linalg.norm(d_vec, axis=(2,1))
            d_norm = np.append(0, d_norm)
            d = np.cumsum(d_norm)
            self._mep_separations = d
            self._mep_length = d[-1]

            xi = d / (d[-1] - d[0])
            self._xi = xi

            # spline the MEP
            self._spline_mep = CubicSpline(x=self._xi, y=self._positions, axis=0)

    def get_samples(
        self,
        temperature: float,
        sampling_mode: str,
        mep_density: float,
        match_density_at_minimum: bool = False,
        n_points_min: int = None
    ) -> dict:
        """Perform thermal sampling around an MEP or a minimum.

        This is the main interface of this class.

        Arguments:
            temperature: temperature in Kelvin for the thermal sampling
            sampling_mode: 'classical' or 'quantum'; select whether effective quantum temperatures
                should be calculated for each mode
            mep_density: linear density of points per unit length (a0) along the MEP
            match_density_at_minimum: if True, the density of points at the minimum is matched 
                with the MEP density without the need to specify `n_points_min`
            n_points_min: number of points to sample at the minimum

        Returns:
            A dictionary of:
                xis: values of xi on the MEP from which the samples originate
                distortions: generated thermal distortions, array of shape (number of distortions, number of atoms, 3)
                labels: integer labels which denote the origin of each distortion, array of shape (number of distortions)
        """
        # check that (only) one of `n_points_min` and `match_density_at_minimum` is specified
        if match_density_at_minimum is False and n_points_min is None:
            raise ValueError('Either `n_points_min` or `match_density_at_minimum` must be specified.')

        # check that we are passing only one of `n_points_min` and `match_density_at_minimum`
        if match_density_at_minimum is True and n_points_min is not None:
            raise ValueError('Only one of `n_points_min` and `match_density_at_minimum` can be used.')

        # set up relevant objects
        idx_nm = self._idx_nm
        xis = []
        distortions = []
        labels = []
        i_label = 0

        # recalculate the mep density to get the density per dimensionless xi 
        scaled_mep_density = mep_density * self._mep_separations[idx_nm[-1]]

        # start sampling from the first minimum geometry
        beta = self._calculate_beta(idx_nm[0], temperature, sampling_mode)
        if match_density_at_minimum:
            n_points_min = self._calculate_minimum_number_of_points(idx=idx_nm[0], density=mep_density, temperature=temperature, beta=beta)
        d_init = self._generate_distortions(
            0,
            np.repeat(self._positions[idx_nm[0]][np.newaxis, ...], n_points_min, axis=0),
            beta
        )

        # if we only have one image, keep all the thermal samples around the image and finish
        if self._number_of_replicas == 1:
            distortions.extend(d_init)
            labels.extend([i_label] * len(d_init))
            xis = None

        # if we have more than one image, then much more stuff needs to be done
        else:
            xi_mep = self._xi
            # filter the thermally sampled minimum to get only points away from the MEP
            d_init = self._process_thermal_samples(d_init, [xi_mep[idx_nm[0]]], operation='filter')
            xis.extend([0.0] * len(d_init))
            distortions.extend(d_init)
            labels.extend([i_label] * len(d_init))
            i_label += 1

            # main loop over the MEP intervals
            for n in range(len(idx_nm[:-1])):

                # inner loop inside each interval
                for i, distribution in enumerate(('cos^2', 'sin^2')):

                    # set up relevant indices: left, right and from where modes are taken in an iteration
                    i_nm_first = idx_nm[n]
                    i_nm_last = idx_nm[n+1]
                    i_nm_modes = idx_nm[n+i]

                    # perform sampling with the chosen parameters
                    xi = random_square_harmonic(distribution, scaled_mep_density, xi_mep[i_nm_first], xi_mep[i_nm_last])
                    beta = self._calculate_beta(i_nm_modes, temperature, sampling_mode)
                    d_mep = self._generate_distortions(i_nm_modes, self._spline_mep(xi), beta)
                    d_proj = self._process_thermal_samples(d_mep, xi, operation='project')  # type: ignore
                    xis.extend(xi)
                    distortions.extend(d_proj)
                    labels.extend([i_label] * len(d_proj))
                    i_label += 1

            # if MEP has modes assigned to the last point, we have to add thermal sampling there, too
            if idx_nm[-1] == (self._number_of_replicas - 1):
                beta = self._calculate_beta(idx_nm[-1], temperature, sampling_mode)
                if match_density_at_minimum:
                    n_points_min = self._calculate_minimum_number_of_points(idx=idx_nm[-1], density=mep_density, temperature=temperature, beta=beta)
                d_fin = self._generate_distortions(
                    idx_nm[-1],
                    np.repeat(self._positions[-1][np.newaxis, ...], n_points_min, axis=0),
                    beta
                )
                d_fin = self._process_thermal_samples(d_fin, [xi_mep[idx_nm[-1]]], operation='filter')
                xis.extend([xi_mep[idx_nm[-1]]] * len(d_fin))
                distortions.extend(d_fin)
                labels.extend([i_label] * len(d_fin))

        samples = {
            'xis': np.array(xis),
            'distortions': np.array(distortions),
            'labels': np.array(labels)
        }

        return samples

    def _process_thermal_samples(self, distortions: Sequence, xi: Sequence, operation: str) -> list:
        """Perform post-generation modification of thermal distortions, either filtering or projecting

        Arguments:
            distortions: sequence of distorted geometries
            xi: sequence of xi values corresponding to the original configuration before distortion
            operation: 'filter' of 'project'
                'filter': only applied at the MEP edges to keep points from the minimum away from MEP
                'project': used along the MEP to project out the direction along the MEP from the distortions

        Returns:
            d_new: modified distortions
        """

        # Note that the `Sequence` type hints abouve might raise some warnings, but it's the sensible thing to do.
        # An issue has been open for... a few years now:
        # https://github.com/numpy/numpy/issues/2776

        # extract stuff from `self` - get coordinates alon MEP and correpsonding tangent vectors
        references = self._spline_mep(xi)
        tangents = self._spline_mep(xi, nu=1)

        # extend `reference` and `tangents` if we only have one for filtering
        # flip derivative if we're at the end of MEP at xi = 1.0
        if operation == 'filter':
            assert len(xi) == 1
            if xi[0] == 1.0:
                tangents *= -1.0
            references = np.repeat(references, len(distortions), axis=0)
            tangents = np.repeat(tangents, len(distortions), axis=0)

        # loop over all distortions, filter or project
        # TODO: possibly broadcast numpy arrays
        d_new = []
        for r, t, d in zip(references, tangents, distortions):
            disp = d - r
            t /= np.linalg.norm(t)
            proj = (disp * t).sum(axis=1).sum(axis=0)
            if operation == 'filter':
                if proj < 0.0:
                    d_new.append(d)
            elif operation == 'project':
                d_proj = d - (t * proj)
                d_new.append(d_proj)
            else:
                raise ValueError('Unknown mode of operation.')

        return d_new

    def _calculate_beta(self, idx: int, temp: float, sampling: str = 'quantum'):
        """Calculate thermal inverse temperature for each mode, classical or quantum.

        Arguments:
            idx: image index, key for the `self._normal_modes` dictionary to select a set of normal modes
            temp: reference temperature in Kelvin
            sampling: either 'classical' or 'quantum'

        Returns:
            beta: array of calculated effective temperatures
        """
        hartree = cs.physical_constants['atomic unit of energy'][0]
        # derive inverse temperatures for the appropriate sampling mode
        if sampling == 'quantum':
            # calculate effective quantum temperatures for each mode (beta* is frequency-dependent)
            if temp > 0.0:
                beta = 1.0 / (cs.k * temp)
                beta *= hartree
                beta = 2 / self._normal_modes[idx].frequencies * np.tanh(beta * self._normal_modes[idx].frequencies / 2)
            elif temp == 0.0:
                # lim beta -> infinity (or as T -> 0) tanh(beta) = 1
                beta = 2 / self._normal_modes[idx].frequencies
            else:
                raise ValueError('Reference classical temperature for quantum sampling must be non-negative.')
        elif sampling == 'classical':
            # stay with just one classical beta, repeat for each frequency to be consistent with the quantum case
            if temp > 0.0:
                beta = 1.0 / (cs.k * temp)
                beta *= hartree
                beta = np.repeat(beta, self._normal_modes[idx].frequencies.size)
            else:
                raise ValueError('Temperature for classical sampling must be positive.')
        else:
            raise ValueError('Unknown mode of sampling.')

        return beta

    def _generate_distortions(self, idx: int, geometries: Iterable, beta: np.ndarray) -> list:
        """Generate multiple distortions, return as a list."""
        distortions = [self._generate_distortion(idx, g, beta) for g in geometries]
        return distortions

    def _generate_distortion(self, idx: int, ref_pos: np.ndarray, beta: np.ndarray) -> np.ndarray:
        """Generate a random distortion from the normal mode thermal distribution.

        In the quantum case, the distribution is a Gaussian just like in the classical case,
        but an effective frequency-dependent temperature
        T* = hbar * omega / (2 * k_B) coth(hbar * omega / (2 * k_B * T))
        is used.

        Arguments:
            idx: image index, key for the `self._normal_modes` dictionary to select a set of normal modes
            ref_pos: reference atomic positions in stationary geometry, shape = (n_atoms, 3)

        Returns:
            Atomic positions of a random thermally distorted structure in nanometers.
        """

        # check `ref_pos` shape and extract number of atoms
        msg = 'Positions need to be a 2D array of shape (n_atoms, 3).'
        assert (ref_pos.ndim == 2) and (ref_pos.shape[1] == 3), msg
        n_atoms = ref_pos.shape[0]

        # reshape everything nicely
        modes = self._normal_modes[idx].modes.reshape(-1, n_atoms, 3)

        # calculate standard deviations of the thermal distributions
        sigmas = 1.0 / (np.sqrt(beta) * self._normal_modes[idx].frequencies)

        # draw random samples of normal coordinates
        q = np.random.normal(loc=0.0, scale=sigmas)

        # prepare distortion vector
        d = (modes * q[:, np.newaxis, np.newaxis]).sum(axis=0)

        # distort the reference structure
        newpositions = ref_pos + d

        return newpositions

    def _calculate_minimum_number_of_points(self, idx: int, density: float, temperature: float, beta: np.ndarray) -> int:
        """Calculate the number of points required to sample the minimum while keeping the density continuous with the MEP sampling density at the minimum.

        Arguments:
            idx: image index, key for the `self._normal_modes` dictionary to select a set of normal modes
            density: linear density of the sampling points in the MEP (per unit length) which we are matching the density of the minimum to
            temperature: temperature of the sampling in K
            beta: array of inverse temperatures, shape = (n_modes)

        Returns:
            n_points: number of points to sample the thermal minimum Gaussian to match the density of the MEP
        """

        # TODO: check that this is correct
        # warn user that this is sort of beta version - should work, but maybe does not
        warnings.warn("""Note that the density matching is a crude implementation that hasn not been tested thoroughly. 
        Use with caution.""")

        # get constants
        amu = cs.physical_constants['atomic mass constant'][0]
        m_e = cs.physical_constants['atomic unit of mass'][0]
        amu_to_me = amu / m_e

        # get key quantities from the class
        xi_min = self._xi[idx]
        masses = np.repeat(self._normal_modes[idx].masses, 3) * amu_to_me
        frequencies = self._normal_modes[idx].frequencies

        # calculate classical beta
        beta_cl = self._calculate_beta(idx=idx, temp=temperature, sampling='classical')[0]

        # get the normal mode vectors in a convenient form
        modes = self._normal_modes[idx].modes
        modes = modes * np.sqrt(masses)
        modes = modes.T

        # get the tangent vector (tau) at the minimum pointing away from the MEP
        tau = self._spline_mep(xi_min, nu=1).reshape(-1)
        if xi_min == 0.0:
            tau = tau * -1.0

        # mass-weight the tangent vector and normalize
        tau = tau * np.sqrt(masses)
        tau = tau / np.linalg.norm(tau)

        # transform tau into the normal mode basis
        tau_omega = modes.T @ tau
        tau_omega = tau_omega / np.linalg.norm(tau_omega)

        # prepare the rest of the tau-oriented basis
        perp = prepare_perp(tau_omega)
        tau_basis_omega = np.append(tau_omega[:, np.newaxis], perp, axis=1)

        # check that the tau-basis has the full rank (i.e. that `tau_omega` is not collinear with other vectors)
        # raise a non-implicit error if this is the case
        # TODO: build the `perp` preparation on random vectors rather than the fixed unit matrix:
        # like that, we can just create a new set if the rank test fails
        if np.linalg.matrix_rank(tau_basis_omega) < tau_basis_omega.shape[0]:
            raise NotImplementedError("""The tau-basis is not full rank. Please use the `n_points_min` argument
            to set the number of points to sample the minimum
            and talk to Krystof to implement random `perp` initialization.""")

        # get the normal-mode-basis diagonal Hessian
        hessian_omega = np.diag((beta / beta_cl) * frequencies**2)

        # transform the Hessian into the tau-oriented basis
        hessian_tau = tau_basis_omega.T @ hessian_omega @ tau_basis_omega

        # calculate the effective squared frequency in the tau direction
        # this uses the formula with determinants I derived - det of the full hessian divided by the det of the minor over the tau subspace
        # we must use the tau-oriented basis, since the minor is basis dependent and only works nicely in this basis
        omega_eff_sq = np.linalg.det(hessian_tau) / np.linalg.det(hessian_tau[1:, 1:])

        # calculate the effective sigma in the tau direction
        sigma_eff = 1.0 / np.sqrt(omega_eff_sq * beta_cl)

        # calculate effective mass of the tau direction
        mu = (tau**2 * masses).sum() / (tau**2).sum()

        # calculate the required number of points
        # alpha is the scaling: alpha * p(tau=0) = density, where p(tau=0) = 1 / (sqrt(2*pi) * sigma_eff)
        alpha = np.sqrt(2 * np.pi) * sigma_eff * mu**-0.5 * density
        n_points = int(alpha)

        return n_points


def random_square_harmonic(
    distribution: str = 'cos^2',
    density: float = 1000.0,
    a: float = 0.0,
    b: float = 1.0
) -> np.ndarray:
    """Generate random samples from a distribution defined by a
    square of harmonic function, either
        cos^2(pi * x / (2 * |b - a|))
    or
        sin^2(pi * x / (2 * |b - a|))
    over the interval [a, b].

    Arguments:
        distribution: 'cos^2' or 'sin^2' to decide whether to sample cosine or sine distribution
        density: required density of samples at the peak of the distribution
        a: lower bound of the independent variable
        b: upper bound of the independent variable

    Returns:
        samples: array of one-dimensional points sampled from the chosen distribution.
    """

    if b <= a:
        raise ValueError('Upper bound must be larger than lower bound.')

    N_samples = density * (b - a)
    N_samples = int(N_samples)

    # points are sampled not from `a` to `b` but from 0 to `b - a` and shifted later by `a`
    points = np.random.uniform([0, 0], [(b - a), 1.0], size=[N_samples, 2])

    if distribution == 'cos^2':
        expected_y = (np.cos(np.pi * points[:, 0] / (2 * (b - a))))**2
    elif distribution == 'sin^2':
        expected_y = (np.sin(np.pi * points[:, 0] / (2 * (b - a))))**2
    else:
        raise ValueError('Unknown distribution type.')

    is_within = points[:, 1] <= expected_y
    samples = points[:, 0][is_within]

    # shift points by a to get correct absolute positions
    samples += a

    return samples


def prepare_perp(vector):
    """Prepare the perpendicular complement (perp) to `vector`
    """
    # scipy orth works only on matrices, not 1D vectors
    assert vector.size > 2

    # initiate perp as identity matrix with the first column taken out (this will be replaced by `vector`)
    perp = np.identity(vector.shape[0])
    perp = perp[:, 1:]

    # project out the direction of `vector` of all remaining columns of `perp`
    proj = np.dot(vector, perp)
    perp -= vector[:, np.newaxis] * proj

    # orthonormalize perp using Scipy SVD
    perp = orth(perp)

    return perp
