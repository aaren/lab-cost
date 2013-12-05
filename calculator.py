import numpy as np

import proc


w = 0.2  # width of tank (m)
Lt = 5


def single_layer_quantities(r=1.01, H=0.25, L=0.25):
    """Calculate quantities of mkp, glycerol needed for
    a single layer RI matched run.

    r - density ratio
    H - depth of fluid
    L - length of lock
    """
    # volumes in L
    v1 = H * w * L * 1000
    v2 = H * w * Lt * 1000

    run = proc.RIMatched(density_ratio=r,
                         sub1='MKP', sub2='Gly',
                         V1=v1, V2=v2)

    quantities = {'MKP': run.mkp.absolute_mass,
                  'Gly': run.gly.absolute_mass}

    return quantities


def density_mix(S, rho1, rho2):
    """Calculate the density of rho0 needed to get
    stratification S given rho1 and rho2, using

        S = (rho1 - rho2) / (rho0 - rho2)
    """
    return (rho1 - rho2) / S + rho2


def two_layer_quantities(S=0.5, h1=0.25, H=0.25, D=0.1, L=0.25, r_min=0.01):
    """Calculate quantities of mkp, glycerol (salt?) needed
    for a two layer run with RI matching in the current and
    lower layer.

    S - stratification parameter
    h1 - depth of lower layer (fraction of H)
    H - depth of fluid
    D - depth of lock
    L - length of lock
    r_min - minimum density difference

    Computes the quantities needed by taking the upper layer as fresh
    water, the lower layer as glycerol solution and the lock fluid as
    a solution of mkp and glycerol.
    """
    # volumes (in L)
    v_lock = w * L * D * 1000
    v_lower = h1 * H * w * Lt * 1000
    # v_upper = (1 - h1) * H * w * Lt * 1000

    # create some aqueous solutions
    gly = proc.AqueousSolution('Gly')
    mkp = proc.AqueousSolution('MKP')

    # range of physical N
    N = np.linspace(1.3333, 1.35)

    # compute necessary lock density over n for given S
    rho_lower = gly.density(N)
    rho_upper = 1  # fresh water
    rho_lock = density_mix(S, rho_lower, rho_upper)

    density_diff_cond = ((rho_lower - rho_upper) > r_min) \
                        & ((rho_lock - rho_lower) > r_min)

    # find the minimum n that satisfies the condition
    n = N[density_diff_cond].min()
    rho_lock = rho_lock[density_diff_cond].min()

    # check that the lock density at this n is not lower than the
    # density of mkp
    mkp = proc.AqueousSolution('MKP')
    if rho_lock > mkp.density(n):
        return 0

    # lower layer mix
    lower_layer = proc.AqueousSolution('Gly', n=n, volume=v_lower)

    # lock consists of two mixtures
    # for given rho_mix, V_mix we can compute what the volumes of
    # these are
    v_lock_gly = v_lock * (mkp.density(n) - rho_lock) \
                        / (mkp.density(n) - gly.density(n))
    v_lock_mkp = v_lock - v_lock_gly

    lock_gly = proc.AqueousSolution('Gly', volume=v_lock_gly, n=n)
    lock_mkp = proc.AqueousSolution('MKP', volume=v_lock_mkp, n=n)

    total_gly = lock_gly.absolute_mass + lower_layer.absolute_mass
    total_mkp = lock_mkp.absolute_mass

    quantities = {'MKP': total_mkp,
                  'Gly': total_gly}

    return quantities


run_matrix = {'single_layer': [{'r': 1.02, 'H': 0.25, 'L': 0.25}] * 7
                              + [{'r': 1.01, 'H': 0.25, 'L': 0.5}] * 5,
              'two_layer':    [{'r_min': 0.02, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.5},
                               {'r_min': 0.02, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.5},
                               {'r_min': 0.02, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.5},
                               {'r_min': 0.02, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.5},
                               {'r_min': 0.01, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.75},
                               {'r_min': 0.01, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.75},
                               {'r_min': 0.01, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.75},
                               {'r_min': 0.01, 'H': 0.25, 'L': 0.5, 'D': 0.25,
                                'h1': 0.25, 'S': 0.75},
                               ]

              }


def compute_totals():
    single_total_gly = sum(single_layer_quantities(**run)['Gly']
                           for run in run_matrix['single_layer'])
    single_total_mkp = sum(single_layer_quantities(**run)['MKP']
                           for run in run_matrix['single_layer'])

    double_total_gly = sum(two_layer_quantities(**run)['Gly']
                           for run in run_matrix['two_layer'])
    double_total_mkp = sum(two_layer_quantities(**run)['MKP']
                           for run in run_matrix['two_layer'])

    total_gly = single_total_gly + double_total_gly
    total_mkp = single_total_mkp + double_total_mkp

    return {'Gly': total_gly, 'MKP': total_mkp}
