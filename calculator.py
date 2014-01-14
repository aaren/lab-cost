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


def density_mix_lock(S, rho1, rho2):
    """Calculate the density of rho0 needed to get
    stratification S given rho1 and rho2, using

        S = (rho1 - rho2) / (rho0 - rho2)
    """
    return (rho1 - rho2) / S + rho2


def density_mix_lower(S, rho0, rho2):
    """Calculate the density of rho1 needed to get
    stratification S given rho0 and rho2, using

        S = (rho1 - rho2) / (rho0 - rho2)
    """
    return (rho0 - rho2) * S + rho2


def two_layer_quantities(S=0.5, h1=0.25, H=0.25, D=0.1, L=0.25, r_min=0.01,
                         v_lock=None, v_lower=None, v_upper=None,
                         t_lock=20, t_lower=20, t_upper=20):
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
    v_lock = v_lock or w * L * D * 1000
    v_lower = v_lower or h1 * H * w * Lt * 1000
    v_upper = v_upper or (1 - h1) * H * w * Lt * 1000

    # create some aqueous solutions
    gly = proc.AqueousSolution('Gly')
    mkp = proc.AqueousSolution('MKP')

    # range of physical refractive index
    N = np.linspace(1.3333, 1.35, 200)

    # compute necessary densities over n for given S
    rho_upper = gly.density(N)
    rho_lock = mkp.density(N)
    rho_lower = density_mix_lower(S, rho_lock, rho_upper)

    density_diff_cond = ((rho_lower - rho_upper) > r_min) \
                        & ((rho_lock - rho_lower) > r_min)

    # find the minimum n that satisfies the condition
    i_min = np.min(np.where(density_diff_cond))
    n = N[i_min]
    rho_lower = rho_lower[i_min]
    rho_upper = rho_upper[i_min]
    rho_lock = rho_lock[i_min]

    # check that the lock density at this n is not lower than the
    # density of mkp
    if rho_lock > mkp.density(n):
        return 0

    # upper layer mix
    upper_layer = proc.AqueousSolution('Gly', n=n, volume=v_upper, temperature=t_upper)
    # lock mix
    lock = proc.AqueousSolution('MKP', n=n, volume=v_lock, temperature=t_lock)

    # lower layer consists of two mixtures. for given rho_mix, V_mix
    # we can compute what the volumes of these are.
    # mass conservation
    v_lower_gly = v_lower * (mkp.density(n) - rho_lower) \
                           / (mkp.density(n) - gly.density(n))
    # volume conservation
    v_lower_mkp = v_lower - v_lower_gly

    lower_gly = proc.AqueousSolution('Gly', volume=v_lower_gly, n=n, temperature=t_lower)
    lower_mkp = proc.AqueousSolution('MKP', volume=v_lower_mkp, n=n, temperature=t_lower)

    total_gly = lower_gly.absolute_mass + upper_layer.absolute_mass
    total_mkp = lower_mkp.absolute_mass + lock.absolute_mass

    print "LOCK"
    print "---------------------"
    print lock.instructions
    print "UPPER LAYER"
    print "---------------------"
    print upper_layer.instructions
    print "LOWER LAYER"
    print "---------------------"
    print lower_gly.instructions
    print lower_mkp.instructions

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
