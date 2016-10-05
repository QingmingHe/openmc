#!/usr/bin/env python
# -*- coding: utf-8 -*-
import openmc


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

    def build_wims69e(self):
        self.group_boundaries = [
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
        self.first_res = 13
        self.last_res = 45


class HomoProblem(object):

    def __init__(self, nuclide, dilution, batches=10, particles=100):

        # Material is composed of H-1 and the object nuclide
        h1 = openmc.Nuclide('H1')
        nuc = openmc.Nuclide(nuclide)
        mat = openmc.Material(material_id=1, name='mat')
        mat.set_density('g/cc', 4.5)
        mat.add_nuclide(h1, dilution / 20.0)
        mat.add_nuclide(nuc, 1.0)
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

        # Set the running parameters
        # energy cutoff??????
        settings_file = openmc.Settings()
        settings_file.run_mode = 'fixed source'
        settings_file.batches = batches
        settings_file.particles = particles
        bounds = [-1, -1, -1, 1, 1, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],
                                        only_fissionable=False)
        watt_dist = openmc.stats.Watt()
        settings_file.source = openmc.source.Source(space=uniform_dist,
                                                    energy=watt_dist)
        settings_file.export_to_xml()


class FullXs(object):

    def __init__(self,
                 temperatures=[294.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0],
                 legendre_order=1,
                 group_structure=GroupStructure('wims69e'),
                 reference_dilution=1e10):
        # Nuclide to be processed should be given by user
        self.nuclide = None

        # Set default settings
        self.temperatures = temperatures
        self.legendre_order = legendre_order
        self.group_structure = group_structure
        self.reference_dilution = reference_dilution

    def generate(self):
        pass


class RItable(object):

    def __init__(self):
        pass


if __name__ == '__main__':
    h = HomoProblem('U238', 28.0)
