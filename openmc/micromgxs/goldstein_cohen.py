# -*- coding: utf-8 -*-
from bisect import bisect_left

_lambda_A = [
    2,
    6,
    11,
    14,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    34,
    35,
    36,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    47,
    48,
    49,
    51,
    52,
    54,
    55,
    57,
    59,
    61,
    62,
    64,
    66,
    69,
    71,
    73,
    75,
    78,
    81,
    84,
    86,
    90,
    93,
    96,
    100,
    103,
    107,
    111,
    116,
    120,
    125,
    131,
    137,
    144,
    152,
    161,
    172,
    188,
    214,
    1000
]

_lambdas = [
    1.00,
    0.99,
    0.98,
    0.97,
    0.96,
    0.95,
    0.94,
    0.93,
    0.91,
    0.90,
    0.89,
    0.87,
    0.85,
    0.84,
    0.82,
    0.81,
    0.79,
    0.78,
    0.76,
    0.75,
    0.74,
    0.72,
    0.71,
    0.70,
    0.69,
    0.68,
    0.67,
    0.66,
    0.65,
    0.64,
    0.63,
    0.62,
    0.61,
    0.60,
    0.59,
    0.58,
    0.57,
    0.56,
    0.55,
    0.54,
    0.53,
    0.52,
    0.51,
    0.50,
    0.49,
    0.48,
    0.47,
    0.46,
    0.45,
    0.44,
    0.43,
    0.42,
    0.41,
    0.40,
    0.39,
    0.38,
    0.37,
    0.36,
    0.35,
    0.34,
    0.33,
    0.32,
    0.31,
    0.30,
    0.29,
    0.28,
    0.27,
    0.26,
    0.25,
    0.24,
    0.23,
    0.22,
    0.21,
    0.20
]


def average_lambda(A):
    i = bisect_left(_lambda_A, A)
    return _lambdas[i]


def _set_u238_lambda():
    import h5py
    import numpy as np
    f = h5py.File('jeff-3.2-wims69e.h5')
    if 'lambda' in f['/U238/resonance']:
        del f['/U238/resonance']['lambda']
    f['/U238/resonance']['lambda'] = np.array([
        1.00000000000E+00, 9.88330500000E-01, 1.00000000000E+00,
        9.71728900000E-01, 9.60144600000E-01, 9.45288900000E-01,
        8.81759600000E-01, 8.98695600000E-01, 5.94728400000E-01,
        4.33874400000E-01, 2.07093700000E-01, 6.38040700000E-02,
        6.45154300000E-02, 4.92320400000E-01, 2.50276300000E-02,
        2.00000000000E-01, 1.00000000000E+00, 1.00000000000E+00,
        2.00000000000E-01, 2.00000000000E-01, 1.00000000000E+00,
        1.00000000000E+00, 1.00000000000E+00, 1.00000000000E+00,
        1.00000000000E+00, 1.00000000000E+00, 1.00000000000E+00,
        1.00000000000E+00, 1.00000000000E+00, 1.00000000000E+00,
        1.00000000000E+00, 1.00000000000E+00, 2.00000000000E-01
    ])
    f.close()

if __name__ == "__main__":
    _set_u238_lambda()
