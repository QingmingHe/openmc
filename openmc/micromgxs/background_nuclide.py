#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import h5py
from math import log10, ceil

K_BOLTZMANN = 8.6173324e-11


class BackgroundNuclide(object):

    def __init__(self, name, potential, A, Z, awr, n_energy_point,
                 start_energy=1e-11, end_energy=20.0,
                 temperatures=[293.6, 600.0, 900.0, 1200.0, 1500.0, 1800.0]):
        self._name = name
        self._potential = potential
        self._start_energy = start_energy
        self._end_energy = end_energy
        self._n_energy_point = n_energy_point
        self._temperatures = temperatures
        self._A = A
        self._Z = Z
        self._awr = awr

        self._temp_strs = [str(int(ceil(temp))) + 'K' for temp in temperatures]
        self._energy = {}
        self._kTs = {}
        self._reactions = {'reaction_002': {}}
        reaction2 = self._reactions['reaction_002']
        for itemp, temp_str in enumerate(self._temp_strs):
            self._energy[temp_str] \
                = np.logspace(log10(start_energy), log10(end_energy),
                              n_energy_point, True)
            self._kTs[temp_str] = K_BOLTZMANN * temperatures[itemp]
            reaction2[temp_str] = np.zeros(n_energy_point) + potential

    def export_to_h5(self, fname=None, fh=None):
        if fh is not None:
            f = fh
        elif fname is not None:
            f = h5py.File(fname)

        nuc_grp = '/' + self._name
        if nuc_grp in f:
            del f[nuc_grp]
        f.create_group(nuc_grp)
        f[nuc_grp].attrs['A'] = self._A
        f[nuc_grp].attrs['Z'] = self._Z
        f[nuc_grp].attrs['atomic_weight_ratio'] = self._awr
        f[nuc_grp].attrs['metastable'] = 0
        erg_grp = nuc_grp + '/energy'
        f.create_group(erg_grp)
        for temp_str in self._temp_strs:
            f[erg_grp][temp_str] = self._energy[temp_str]
        kTs_grp = nuc_grp + '/kTs'
        f.create_group(kTs_grp)
        for temp_str in self._kTs:
            f[kTs_grp][temp_str] = self._kTs[temp_str]
        reactions_grp = nuc_grp + '/reactions'
        f.create_group(reactions_grp)
        for reaction in self._reactions:
            reaction_grp = reactions_grp + '/' + reaction
            f.create_group(reaction_grp)
            f[reaction_grp].attrs['Q_value'] = 0.0
            f[reaction_grp].attrs['center_of_mass'] = 1
            f[reaction_grp].attrs['label'] = np.string_('(n,elastic)')
            f[reaction_grp].attrs['mt'] = 2
            for temp_str in self._temp_strs:
                rea_temp_grp = reaction_grp + '/' + temp_str
                f.create_group(rea_temp_grp)
                f[rea_temp_grp]['xs'] \
                    = self._reactions[reaction][temp_str]
                f[rea_temp_grp]['xs'].attrs['threshold_idx'] = 1

        if fname is not None:
            f.close()

if __name__ == '__main__':
    nuc = BackgroundNuclide('H1b', 20.0, 1, 1, 0.999167, 100)
    nuc.export_to_h5('background.h5')
