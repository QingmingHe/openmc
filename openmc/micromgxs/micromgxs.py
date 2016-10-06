#!/usr/bin/env python
# -*- coding: utf-8 -*-
import openmc
import os
import numpy as np
from glob import glob


class MicroMgXs(object):

    def __init__(self):
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
        self._first_res = 13
        self._last_res = 45

    @property
    def group_boundaries(self):
        return self._group_boundaries

    @property
    def ng(self):
        return len(self._group_boundaries) - 1


def export_homo_problem_xml(nuclide, dilution, temperature,
                            background_nuclide, fission_nuclide=None,
                            fisnuc_refdil=None):
    # Material is composed of H-1 and the object nuclide
    h1 = openmc.Nuclide(background_nuclide)
    nuc = openmc.Nuclide(nuclide)
    mat = openmc.Material(material_id=1, name='mat')
    mat.set_density('g/cc', 4.5)
    mat.add_nuclide(h1, dilution / 20.0)
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
        export_homo_problem_xml(self._nuclide, self._reference_dilution,
                                temperature, self._background_nuclide)

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

        energy_tally = openmc.Tally()
        grp_bnds = [val / 1e6 for val in
                    sorted(self._group_structure.group_boundaries)]
        energy_filter = openmc.EnergyFilter(grp_bnds)
        energy_out_filter = openmc.EnergyoutFilter(grp_bnds)
        cell_filter = openmc.CellFilter((1, ))

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
        flux_tally = openmc.Tally()
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

        flux_tally.estimator = 'tracklength'
        flux_tally.filters = [energy_filter, cell_filter]
        flux_tally.scores = ['flux']

        tallies = openmc.Tallies()
        tallies.append(energy_out_tally)
        tallies.append(flux_tally)
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
        self._nu_fission_matrix = np.zeros((ng, ng))
        self._flux_fix = np.zeros((n_temp, ng))
        self._flux_ein = np.zeros(ng)
        self._nu_scatter = np.zeros((n_temp, ng_in, ng_out,
                                     self._legendre_order + 1))
        self._chi = np.zeros(ng)

        # for temperature in self._temperatures:
        for itemp, temperature in enumerate(self._temperatures):
            # Export fixed source input files
            self._export_fs_xml(temperature)

            # Run OpenMC
            openmc.run()

            # Load the tally data from statepoint
            statepoint = glob(os.path.join(
                os.getcwd(), "statepoint.%s.*" % (self._batches)))[0]
            self._load_fix_statepoint(statepoint, itemp)

            # Export eigenvalue input files
            self._export_eig_xml()

            # Run OpenMC
            # openmc.run()

            # Load the tally data from statepoint
            # statepoint = glob(os.path.join(
            #     os.getcwd(), "statepoint.%s.*" % (self._batches)))[0]

    def _load_fix_statepoint(self, statepoint, itemp):
        sp = openmc.StatePoint(statepoint)
        ng = self._group_structure.ng

        # Get flux
        self._flux_fix[itemp, :] = sp.get_tally(scores=['flux'])\
                                     .sum[:, 0, 0][::-1]

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
        pass


if __name__ == '__main__':
    f = FullXs()
    f.nuclide = 'U238'
    f.reference_dilution = 28.0
    f.build_library()
