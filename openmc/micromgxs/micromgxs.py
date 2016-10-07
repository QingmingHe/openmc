#!/usr/bin/env python
# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import openmc
import os
import numpy as np
from glob import glob
import h5py
from math import pi
import re

RESONANCE_FISSION_AUTO = 1
RESONANCE_FISSION_USER = 2


def _get_potentials_from_endf(endf_path, potentials_fname):
    files = glob(os.path.join(endf_path, "*"))
    potentials = {}
    names = []
    for afile in sorted(files):
        name = _endf_fname_to_name(afile)
        potential = _get_potential_from_endf(afile)
        potentials[name] = potential
        names.append(name)

    with open(potentials_fname, 'w') as f:
        for name in names:
            f.write('\"%s\": %f,\n' % (name, potentials[name]))


def _endf_fname_to_name(fname):
    basename = os.path.basename(fname)
    names = re.findall("n-\d+-([a-zA-Z]{1,2})-(\d+)([mM])?", basename)[0]
    name = names[0] + str(int(names[1]))
    if len(names[2]) != 0:
        name += "_m1"
    return name


def _float_fortran(string):
    x, y = re.split('[+-]', string)
    if '-' in string:
        return float(x) * 10 ** -float(y)
    else:
        return float(x) * 10 ** float(y)


def _get_potential_from_endf(fname):
    with open(fname) as f:
        for aline in f:
            if aline[71:75] == '2151':
                f.next()
                f.next()
                aline = f.next()
                if aline[12:13] == ' ':
                    return 0.0
                else:
                    a = _float_fortran(aline.strip().split()[1])
                    return 4.0 * pi * a ** 2


def _has_resfis(A, Z):
    has_resfis = False
    if A == 92 and (Z == 233 or Z == 235):
        has_resfis = True
    elif A == 94 and (Z == 239 or Z == 241):
        has_resfis = True
    return has_resfis


def _get_A_Z_awr(cross_sections, materials):
    tree = ET.parse(cross_sections)
    root = tree.getroot()
    direc = os.path.dirname(cross_sections)
    for child in root:
        if materials == child.attrib['materials']:
            path = os.path.join(direc, child.attrib['path'])
            f = h5py.File(path, 'r')
            A = f[materials].attrs['A']
            Z = f[materials].attrs['Z']
            awr = f[materials].attrs['atomic_weight_ratio']
            f.close()
            return A, Z, awr
    raise Exception('%s cannot be found in %s!' % (materials, cross_sections))


def _condense_scatter(x):
    ng = len(x)
    for ig0, val in enumerate(x):
        if val != 0.0:
            break
    for ig1, val in enumerate(x[::-1]):
        if val != 0.0:
            break
    ig1 = ng - ig1
    return ig0, ig1, x[ig0:ig1]


class MicroMgXsLibrary(object):

    def __init__(self):
        pass

    def export_to_h5(self, fname):
        # group structure
        pass


class MicroMgXsNuclide(object):

    def __init__(self):
        pass

    def export_to_h5(self, fname):
        # A, Z, awr
        pass


class GroupStructure(object):

    def __init__(self, name=None, group_boundaries=None, first_res=None,
                 last_res=None):
        if name is not None:
            if name == 'wims69e':
                self.build_wims69e()
            else:
                raise Exception('group structure name %s is not supported!' %
                                (name))
        else:
            if group_boundaries is None:
                raise Exception('group_boundaries should be given!')
            if first_res is None:
                raise Exception('first_res should be given!')
            if last_res is None:
                raise Exception('last_res should be given!')
            self._group_boundaries = group_boundaries
            self._first_res = first_res
            self._last_res = last_res

    def build_wims69e(self):
        self._group_boundaries = [
            1.00000E+07, 6.06550E+06, 3.67900E+06, 2.23100E+06, 1.35300E+06,
            8.21000E+05, 5.00000E+05, 3.02500E+05, 1.83000E+05, 1.11000E+05,
            6.73400E+04, 4.08500E+04, 2.47800E+04, 1.50300E+04, 9.11800E+03,
            5.53000E+03, 3.51910E+03, 2.23945E+03, 1.42510E+03, 9.06899E+02,
            3.67263E+02, 1.48729E+02, 7.55014E+01, 4.80520E+01, 2.77000E+01,
            1.59680E+01, 9.87700E+00, 4.00000E+00, 3.30000E+00, 2.60000E+00,
            2.10000E+00, 1.50000E+00, 1.30000E+00, 1.15000E+00, 1.12300E+00,
            1.09700E+00, 1.07100E+00, 1.04500E+00, 1.02000E+00, 9.96000E-01,
            9.72000E-01, 9.50000E-01, 9.10000E-01, 8.50000E-01, 7.80000E-01,
            6.25000E-01, 5.00000E-01, 4.00000E-01, 3.50000E-01, 3.20000E-01,
            3.00000E-01, 2.80000E-01, 2.50000E-01, 2.20000E-01, 1.80000E-01,
            1.40000E-01, 1.00000E-01, 8.00000E-02, 6.70000E-02, 5.80000E-02,
            5.00000E-02, 4.20000E-02, 3.50000E-02, 3.00000E-02, 2.50000E-02,
            2.00000E-02, 1.50000E-02, 1.00000E-02, 5.00000E-03, 1.00000E-05]
        self._first_res = 12
        self._last_res = 45

    @property
    def res_group_bnds(self):
        return self._group_boundaries[self._first_res:self._last_res+1]

    @property
    def group_boundaries(self):
        return self._group_boundaries

    @property
    def first_res(self):
        return self._first_res

    @property
    def last_res(self):
        return self._last_res

    @property
    def ng(self):
        return len(self._group_boundaries) - 1

    @property
    def n_res(self):
        return self._last_res - self._first_res


def export_homo_problem_xml(nuclide, dilution, temperature,
                            background_nuclide, fission_nuclide=None,
                            fisnuc_refdil=None):
    # Material is composed of H-1 and the object nuclide
    h1 = openmc.Nuclide(background_nuclide)
    nuc = openmc.Nuclide(nuclide)
    mat = openmc.Material(material_id=1, name='mat')
    mat.set_density('atom/b-cm', 0.069335)
    mat.add_nuclide(h1, dilution / 20.478)
    mat.add_nuclide(nuc, 1.0)
    if fission_nuclide is not None:
        if fission_nuclide != nuclide:
            if fisnuc_refdil is None:
                raise Exception(
                    'fisnuc_refdil should be given')
            fisnuc = openmc.Nuclide(fission_nuclide)
            mat.add_nuclide(fisnuc, dilution / fisnuc_refdil)
    materials_file = openmc.Materials([mat])
    materials_file.export_to_xml()

    # Cell is box with reflective boundary
    x1 = openmc.XPlane(surface_id=1, x0=-1)
    x2 = openmc.XPlane(surface_id=2, x0=1)
    y1 = openmc.YPlane(surface_id=3, y0=-1)
    y2 = openmc.YPlane(surface_id=4, y0=1)
    z1 = openmc.ZPlane(surface_id=5, z0=-1)
    z2 = openmc.ZPlane(surface_id=6, z0=1)
    for surface in [x1, x2, y1, y2, z1, z2]:
        surface.boundary_type = 'reflective'
    box = openmc.Cell(cell_id=1, name='box')
    box_region = +x1 & -x2 & +y1 & -y2 & +z1 & -z2
    box.region = box_region
    box.fill = mat
    root = openmc.Universe(universe_id=0, name='root universe')
    root.add_cell(box)
    geometry = openmc.Geometry(root)
    geometry.export_to_xml()

    return mat, geometry


class FullXs(object):

    def __init__(self):
        # Nuclide to be processed should be given by user
        self._nuclide = None

        # Set default settings
        self._temperatures = [294.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0]
        self._legendre_order = 1
        self._group_structure = GroupStructure('wims69e')
        self._reference_dilution = 1e10
        self._background_nuclide = 'H1'
        self._fission_nuclide = 'U235'
        self._fisnuc_refdil = 800.0
        self._batches = 10
        self._particles = 100
        self._inactive = 5

    def _export_fs_xml(self, temperature):
        # Export geometry and materials of homogeneous problem
        self._material, self._geometry = export_homo_problem_xml(
            self._nuclide, self._reference_dilution, temperature,
            self._background_nuclide)

        # Calculate number density of object nuclide
        sum_dens = 0.0
        nuclides = self._material.get_nuclide_densities()
        for nuc in nuclides:
            sum_dens += nuclides[nuc][1]
        self._nuclide_density \
            = self._material._density * nuclides[self._nuclide][1] / sum_dens

        # Set the running parameters
        # energy cutoff??????
        settings_file = openmc.Settings()
        settings_file.run_mode = 'fixed source'
        settings_file.batches = self._batches
        settings_file.particles = self._particles
        bounds = [-1, -1, -1, 1, 1, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],
                                        only_fissionable=False)
        watt_dist = openmc.stats.Watt()
        settings_file.source = openmc.source.Source(space=uniform_dist,
                                                    energy=watt_dist)
        settings_file.temperature = {'default': temperature}
        settings_file.export_to_xml()

        # Create tallies
        tallies = openmc.Tallies()

        grp_bnds = [val / 1e6 for val in
                    sorted(self._group_structure.group_boundaries)]
        energy_filter = openmc.EnergyFilter(grp_bnds)
        energy_out_filter = openmc.EnergyoutFilter(grp_bnds)
        cell_filter = openmc.CellFilter((1, ))

        energy_tally = openmc.Tally()
        energy_tally.estimator = 'tracklength'
        energy_tally.filters = [energy_filter, cell_filter]
        energy_tally.nuclides = [self._nuclide]
        energy_tally.scores = ['total', 'fission', 'nu-fission']
        tallies.append(energy_tally)

        for i in range(self._legendre_order + 1):
            energy_out_tally = openmc.Tally()
            energy_out_tally.estimator = 'analog'
            energy_out_tally.filters = [energy_filter, energy_out_filter,
                                        cell_filter]
            energy_out_tally.nuclides = [self._nuclide]
            energy_out_tally.scores = ['nu-scatter-%s' % (i)]
            tallies.append(energy_out_tally)

        flux_tally = openmc.Tally()
        flux_tally.estimator = 'tracklength'
        flux_tally.filters = [energy_filter, cell_filter]
        flux_tally.scores = ['flux']
        tallies.append(flux_tally)

        tallies.export_to_xml()

    def _export_eig_xml(self):
        # Export geometry and materials of homogeneous problem
        export_homo_problem_xml(self._nuclide, self._reference_dilution,
                                self._temperatures[0],
                                self._background_nuclide,
                                fission_nuclide=self._fission_nuclide,
                                fisnuc_refdil=self._fisnuc_refdil)

        # Set the running parameters
        settings_file = openmc.Settings()
        settings_file.run_mode = 'eigenvalue'
        settings_file.batches = self._batches
        settings_file.inactive = self._inactive
        settings_file.particles = self._particles
        bounds = [-1, -1, -1, 1, 1, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],
                                        only_fissionable=False)
        settings_file.source = openmc.source.Source(space=uniform_dist)
        settings_file.temperature = {'default': self._temperatures[0]}
        settings_file.export_to_xml()

        # Create tallies
        energy_out_tally = openmc.Tally()
        grp_bnds = [val / 1e6 for val in
                    sorted(self._group_structure.group_boundaries)]
        energy_filter = openmc.EnergyFilter(grp_bnds)
        energy_out_filter = openmc.EnergyoutFilter(grp_bnds)
        cell_filter = openmc.CellFilter((1, ))
        energy_out_tally.estimator = 'analog'
        energy_out_tally.filters = [energy_filter, energy_out_filter,
                                    cell_filter]
        energy_out_tally.nuclides = [self._nuclide]
        energy_out_tally.scores = ['nu-fission']

        tallies = openmc.Tallies()
        tallies.append(energy_out_tally)
        tallies.export_to_xml()

    def build_library(self):
        # Initialize multi-group cross sections (full xs part)
        ng = self._group_structure.ng
        ng_in = ng
        ng_out = ng
        n_temp = len(self._temperatures)
        self._fissionable = False
        self._total = np.zeros((n_temp, ng))
        self._fission = np.zeros((n_temp, ng))
        self._nu_fission = np.zeros((n_temp, ng))
        self._nu_fission_matrix = np.zeros((ng_in, ng_out))
        self._flux_fix = np.zeros((n_temp, ng))
        self._nu_scatter = np.zeros((n_temp, ng_in, ng_out,
                                     self._legendre_order + 1))
        self._chi = np.zeros(ng)

        for itemp, temperature in enumerate(self._temperatures):
            # Export fixed source input files
            self._export_fs_xml(temperature)

            # Run OpenMC
            openmc.run()

            # Load the tally data from statepoint
            statepoint = glob(os.path.join(
                os.getcwd(), "statepoint.%s.*" % (self._batches)))[0]
            self._load_fix_statepoint(statepoint, itemp)

        if self._fissionable:
            # Process chi for fissionable nuclide
            # Export eigenvalue input files
            self._export_eig_xml()

            # Run OpenMC
            openmc.run()

            # Load the tally data from statepoint
            statepoint = glob(os.path.join(
                os.getcwd(), "statepoint.%s.*" % (self._batches)))[0]
            self._load_eig_statepoint(statepoint)

    def export_to_h5(self, fname):
        f = h5py.File(fname)
        root_group = '/' + self._nuclide
        ng = self._group_structure.ng
        n_temp = len(self._temperatures)

        # Remove existing full xs data
        group = root_group + '/full_xs'
        if group in f:
            del f[group]
        f.create_group(group)

        # Whether fissionable
        if self._fissionable:
            f[group].attrs['fissionable'] = 1
        else:
            f[group].attrs['fissionable'] = 0

        # Legendre order
        f[group].attrs['legendre_order'] = self._legendre_order

        # Temperatures
        f[group + '/temperatures'] = self._temperatures

        for itemp in range(n_temp):
            group = '%s/full_xs/%sK' % (
                root_group, int(self._temperatures[itemp]))

            # Total cross sections
            f[group + '/total'] = self._total[itemp, :]

            # Scattering matrix
            for il in range(self._legendre_order + 1):
                for ig_out in range(ng):
                    ig0, ig1, scatter \
                        = _condense_scatter(
                            self._nu_scatter[itemp, :, ig_out, il])
                    dset = group + '/scatter/%s/%s' % (il, ig_out)
                    f[dset] = scatter
                    f[dset].attrs['ig0'] = ig0
                    f[dset].attrs['ig1'] = ig1

            if self._fissionable:
                # Fission chi
                f[group + '/chi'] = self._chi

                # Fission cross sections
                f[group + '/fission'] = self._fission[itemp, :]

                # Nu fission cross sections
                f[group + '/nu_fission'] = self._nu_fission[itemp, :]

        f.close()

    def _load_fix_statepoint(self, statepoint, itemp):
        sp = openmc.StatePoint(statepoint)
        ng = self._group_structure.ng

        # Get flux
        self._flux_fix[itemp, :] \
            = sp.get_tally(scores=['flux']) .sum[:, 0, 0][::-1] \
            * self._nuclide_density

        # Get total xs
        self._total[itemp, :] = sp.get_tally(scores=['total'])\
                                  .sum[:, 0, 0][::-1]
        self._total[itemp, :] /= self._flux_fix[itemp, :]

        # Get fission xs
        self._fission[itemp, :] = sp.get_tally(scores=['fission'])\
                                    .sum[:, 0, 0][::-1]
        self._fission[itemp, :] /= self._flux_fix[itemp, :]

        # Get nu fission xs
        self._nu_fission[itemp, :] = sp.get_tally(
            scores=['nu-fission']).sum[:, 0, 0][::-1]
        self._nu_fission[itemp, :] /= self._flux_fix[itemp, :]

        # Get scattering matrix
        for i in range(self._legendre_order + 1):
            self._nu_scatter[itemp, :, :, i] \
                = sp.get_tally(scores=['nu-scatter-%s' % i])\
                    .sum[:, 0, 0][::-1].reshape(ng, ng)
            for ig in range(ng):
                self._nu_scatter[itemp, ig, :, i] /= self._flux_fix[itemp, ig]

        # Determine whether is fissionable nuclide
        if sum(self._fission[itemp, :]) != 0.0:
            self._fissionable = True

    def _load_eig_statepoint(self, statepoint):
        sp = openmc.StatePoint(statepoint)
        ng = self._group_structure.ng

        # Get fission matrix
        self._nu_fission_matrix[:, :] = sp.get_tally(scores=['nu-fission'])\
                                          .sum[:, 0, 0][::-1].reshape(ng, ng)

        # Calculate fission chi
        self._chi[:] = self._nu_fission_matrix.sum(axis=0) \
            / self._nu_fission_matrix.sum()

    @property
    def nuclide(self):
        return self._nuclide

    @nuclide.setter
    def nuclide(self, nuclide):
        self._nuclide = nuclide

    @property
    def temperatures(self):
        return self._temperatures

    @temperatures.setter
    def temperatures(self, temperatures):
        self._temperatures = temperatures

    @property
    def legendre_order(self):
        return self._legendre_order

    @legendre_order.setter
    def legendre_order(self, legendre_order):
        self._legendre_order = legendre_order

    @property
    def group_structure(self):
        return self._group_structure

    @group_structure.setter
    def group_structure(self, group_structure):
        if not isinstance(group_structure, GroupStructure):
            raise Exception(
                'group_structure should be instance of GroupStructure')
        self._group_structure = group_structure

    @property
    def reference_dilution(self):
        return self._reference_dilution

    @reference_dilution.setter
    def reference_dilution(self, reference_dilution):
        self._reference_dilution = reference_dilution

    @property
    def background_nuclide(self):
        return self._background_nuclide

    @background_nuclide.setter
    def background_nuclide(self, background_nuclide):
        self._background_nuclide = background_nuclide

    @property
    def batches(self):
        return self._batches

    @batches.setter
    def batches(self, batches):
        self._batches = batches

    @property
    def particles(self):
        return self._particles

    @particles.setter
    def particles(self, particles):
        self._particles = particles

    @property
    def inactive(self):
        return self._inactive

    @inactive.setter
    def inactive(self, inactive):
        self._inactive = inactive

    @property
    def fission_nuclide(self):
        return self._fission_nuclide

    @fission_nuclide.setter
    def fission_nuclide(self, fission_nuclide):
        self._fission_nuclide = fission_nuclide

    @property
    def fisnuc_refdil(self):
        return self._fisnuc_refdil

    @fisnuc_refdil.setter
    def fisnuc_refdil(self, fisnuc_refdil):
        self._fisnuc_refdil = fisnuc_refdil


class RItable(object):

    def __init__(self):
        self._nuclide = None

        # Set default settings
        self._dilutions = [5.0, 1e1, 15.0, 20.0, 25.0, 28.0, 30.0, 35.0, 40.0,
                           45.0, 50.0, 52.0, 60.0, 70.0, 80.0, 1e2, 2e2, 4e2,
                           6e2, 1e3, 1.2e3, 1e4, 1e10]
        self._temperatures = [294.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0]
        self._group_structure = GroupStructure('wims69e')
        self._background_nuclide = 'H1'
        self._batches = 10
        self._particles = 100
        self._resfis_method = RESONANCE_FISSION_AUTO
        self._has_resfis = False

    def _export_xml(self, temperature, dilution):
        # Export geometry and materials of homogeneous problem
        self._material, self._geometry = export_homo_problem_xml(
            self._nuclide, dilution, temperature, self._background_nuclide)

        # Calculate number density of object nuclide
        sum_dens = 0.0
        nuclides = self._material.get_nuclide_densities()
        for nuc in nuclides:
            sum_dens += nuclides[nuc][1]
        self._nuclide_density \
            = self._material._density * nuclides[self._nuclide][1] / sum_dens

        # Set the running parameters
        # energy cutoff??????
        settings_file = openmc.Settings()
        settings_file.run_mode = 'fixed source'
        settings_file.batches = self._batches
        settings_file.particles = self._particles
        bounds = [-1, -1, -1, 1, 1, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],
                                        only_fissionable=False)
        watt_dist = openmc.stats.Watt()
        settings_file.source = openmc.source.Source(space=uniform_dist,
                                                    energy=watt_dist)
        settings_file.temperature = {'default': temperature}
        settings_file.export_to_xml()

        # Create tallies
        tallies = openmc.Tallies()

        grp_bnds = [val / 1e6 for val in
                    sorted(self._group_structure.res_group_bnds)]
        energy_filter = openmc.EnergyFilter(grp_bnds)
        cell_filter = openmc.CellFilter((1, ))

        flux_tally = openmc.Tally()
        flux_tally.estimator = 'tracklength'
        flux_tally.filters = [energy_filter, cell_filter]
        flux_tally.scores = ['flux']
        tallies.append(flux_tally)

        reaction_tally = openmc.Tally()
        reaction_tally.estimator = 'tracklength'
        reaction_tally.filters = [energy_filter, cell_filter]
        reaction_tally.nuclides = [self._nuclide]
        reaction_tally.scores = ['absorption', 'scatter']
        if self._has_resfis:
            reaction_tally.scores.append('nu-fission')
        tallies.append(reaction_tally)

        background_tally = openmc.Tally()
        background_tally.estimator = 'tracklength'
        background_tally.filters = [energy_filter, cell_filter]
        background_tally.nuclides = [self._background_nuclide]
        background_tally.scores = ['scatter']
        tallies.append(background_tally)

        tallies.export_to_xml()

    def build_library(self):
        # Whether has resonance fission
        if self._resfis_method == RESONANCE_FISSION_AUTO:
            cross_sections = os.getenv('OPENMC_CROSS_SECTIONS')
            if cross_sections is None:
                raise Exception('OPENMC_CROSS_SECTIONS env var should be set!')
            self._A, self._Z, self._awr \
                = _get_A_Z_awr(cross_sections, self._nuclide)
            self._has_resfis = _has_resfis(self._A, self._Z)

        # Initialize multi-group cross sections (resonance xs table part)
        n_res = self._group_structure.n_res
        n_temp = len(self._temperatures)
        n_dilution = len(self._dilutions)
        self._flux = np.zeros((n_temp, n_res, n_dilution))
        self._absorption = np.zeros((n_temp, n_res, n_dilution))
        self._nu_scatter = np.zeros((n_temp, n_res, n_dilution))
        self._dilus_group = np.zeros((n_temp, n_res, n_dilution))
        if self._has_resfis:
            self._nu_fission = np.zeros((n_temp, n_res, n_dilution))

        for itemp, temperature, in enumerate(self._temperatures):
            for idil, dilution in enumerate(self._dilutions):
                # Export input files
                self._export_xml(temperature, dilution)

                # Run OpenMC
                openmc.run()

                # Load the tally data from statepoint
                statepoint = glob(os.path.join(
                    os.getcwd(), "statepoint.%s.*" % (self._batches)))[0]
                self._load_statepoint(statepoint, itemp, idil)

    def _load_statepoint(self, statepoint, itemp, idil):
        sp = openmc.StatePoint(statepoint)
        dens1 = self._nuclide_density

        # Get flux
        self._flux[itemp, :, idil] \
            = sp.get_tally(scores=['flux']).sum[:, 0, 0][::-1] * dens1

        # Get absorption xs
        self._absorption[itemp, :, idil] \
            = sp.get_tally(scores=['absorption'], nuclides=[self._nuclide])\
                .sum[:, 0, 0][::-1]
        self._absorption[itemp, :, idil] /= self._flux[itemp, :, idil]

        # Get scatter xs
        self._nu_scatter[itemp, :, idil] \
            = sp.get_tally(scores=['scatter'], nuclides=[self._nuclide])\
                .sum[:, 0, 0][::-1]
        self._nu_scatter[itemp, :, idil] /= self._flux[itemp, :, idil]

        # Get nu fission xs
        if self._has_resfis:
            self._nu_fission[itemp, :, idil] \
                = sp.get_tally(scores=['nu-fission'],
                               nuclides=[self._nuclide]).sum[:, 0, 0][::-1]
            self._nu_fission[itemp, :, idil] \
                /= self._flux[itemp, :, idil]

        # Get dilution xs
        self._dilus_group[itemp, :, idil] \
            = sp.get_tally(scores=['scatter'],
                           nuclides=[self._background_nuclide])\
                .sum[:, 0, 0][::-1]
        self._dilus_group[itemp, :, idil] /= self._flux[itemp, :, idil]

    def export_to_h5(self, fname):
        f = h5py.File(fname)
        root_group = '/' + self._nuclide

        # Create resonance group
        group = root_group + '/resonance'
        if group in f:
            del f[group]
        f.create_group(group)

        # Whether has resonance fission
        if self._has_resfis:
            f[group].attrs['resfis'] = 1
        else:
            f[group].attrs['resfis'] = 0

        # Resonance energy range
        f[group].attrs['first_res'] = self._group_structure.first_res
        f[group].attrs['last_res'] = self._group_structure.last_res

        # Temperatures
        f[group + '/temperatures'] = self._temperatures

        # Resonance cross sections
        f[group + '/absorption'] = self._absorption
        f[group + '/scatter'] = self._nu_scatter
        if self._has_resfis:
            f[group + '/nu_fission'] = self._nu_fission
        f[group + '/dilutions'] = self._dilus_group

        f.close()

    @property
    def has_resfis(self):
        return self._has_resfis

    @has_resfis.setter
    def has_resfis(self, has_resfis):
        self._has_resfis = has_resfis

    @property
    def resfis_method(self):
        return self._resfis_method

    @resfis_method.setter
    def resfis_method(self, resfis_method):
        self._resfis_method = resfis_method

    @property
    def particles(self):
        return self._particles

    @particles.setter
    def particles(self, particles):
        self._particles = particles

    @property
    def batches(self):
        return self._batches

    @batches.setter
    def batches(self, batches):
        self._batches = batches

    @property
    def background_nuclide(self):
        return self._background_nuclide

    @background_nuclide.setter
    def background_nuclide(self, background_nuclide):
        self._background_nuclide = background_nuclide

    @property
    def group_structure(self):
        return self._group_structure

    @group_structure.setter
    def group_structure(self, group_structure):
        self._group_structure = group_structure

    @property
    def temperatures(self):
        return self._temperatures

    @temperatures.setter
    def temperatures(self, temperatures):
        self._temperatures = temperatures

    @property
    def nuclide(self):
        return self._nuclide

    @nuclide.setter
    def nuclide(self, nuclide):
        self._nuclide = nuclide

    @property
    def dilutions(self):
        return self._dilutions

    @dilutions.setter
    def dilutions(self, dilutions):
        self._dilutions = dilutions

if __name__ == '__main__':
    # f = FullXs()
    # f.nuclide = 'U238'
    # f.reference_dilution = 28.0
    # f.build_library()
    # f.export_to_h5(f.nuclide + '.h5')
    # r = RItable()
    # r.nuclide = "U238"
    # r.build_library()
    # r.export_to_h5(r.nuclide + '.h5')

    get_potentials_from_endf("/home/qingming/library/jeff-3.2/", "potentials")
