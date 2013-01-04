from __future__ import division

from sys import argv

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def get_data():
    """Grab crc data from file."""
    fname = 'crc-data'
    d = np.loadtxt(fname)  # ignores lines beginning '#'
    # specify substances and the line index at which they occur
    # (ignoring comments)
    substances = [('Gly', 0, 15),
                  ('NaCl', 15, 34),
                  ('MKP', 34, 46),
                  ('DKP', 46, 59),
                  ('KCl', 59, 77),
                  ('MNP', 77, 107),
                  ('DNP', 107, 118),
                  ('MgCl', 118, 139),
                  ('LiCl', 139, 159),
                  ('EDTA', 159, 170),
                  ('KBr', 170, -1)]
    solubilities = [('Gly', 9999),      # units g substance / L water. (wikipedia)
                    ('NaCl', 359),
                    ('MKP', 220),
                    ('DKP', 150),
                    ('KCl', 281),
                    ('MNP', 599),
                    ('DNP', 77),
                    ('MgCl', 543),
                    ('LiCl', 000),      # FIXME: correct value
                    ('EDTA', 000),       # FIXME: correct value
                    ('KBr', 000)]       # FIXME: correct value
    fields = [('wt.', 0), ('density', 3), ('n', 4), ('viscosity', 6)]
    data = {s: {f: d[a:b, c] for f, c in fields} for s, a, b in substances}
    for sub, sol in solubilities:
        data[sub]['solubility'] = sol
    return data


def max_density(substance, d=None):
    """Calculate the maximum density that a substance can be
    mixed to in water.
    """
    if not d:
        d = get_data()
    m, c = calc_coefficients(d, substance, 'wt.', 'density')

    sol = d[substance]['solubility']
    wt = (sol / (sol + 1000)) * 100
    max_density = m * wt + c
    return max_density


def calc_coefficients(data, substance, x, y):
    """Assuming a straight line fit, calculate (m,c) for specific
    values of x, y (density, n, etc.) for a specific substance.
    """
    m, c = np.polyfit(data[substance][x], data[substance][y], 1)
    return m, c


def n_sensitivity(substance, n, volume=None, dn=0.0001):
    """A prototype experimental n-matching procedure is to mix the
    Glycerol phase to the required density and then measure the
    refractive index, n_gly. The mkp phase is then mixed such that
    n_mkp == n_gly and the density recorded.

    What we want to know is how accurately we can set n_mkp; i.e. for
    the total volume of the lock, what is the quantity of mkp that
    corresponds to dn = 0.0001? If we can't set n to within dn, we are
    a bit stuck.

    Essentially, we can *measure* n to within dn but if we can't
    *control* it to the same precision then this method won't work.

    This function calculates the required mass precision for a given
    substance at a given n and volume, with dn=0.0001 by default.

    Returns: dm, the mass precision.
    """
    # volume for the substance specified
    volumes = {'Gly': 0.2 * 0.25 * 5.5, 'MKP': 0.2 * 0.25 * 0.25}
    if not volume:
        volume = volumes[substance]

    d = get_data()
    rho = density(n, substance, d)
    # total mass
    M = volume * rho * 1000
    # variation of n with % wt.
    m, c = calc_coefficients(d, substance, 'wt.', 'n')
    dwt = dn / m
    dm = M * dwt / 100

    print rho
    print volume
    print M
    print dwt
    print dm

    return dm


def cost(ratio):
    """Takes desired density difference or ratio as input and
    calculates the cost of a run based on unit costs of glycerol
    and mkp.
    Outputs: cost of run, quantities of gly/mkp, required n.
    """
    # determine coefficients
    d = get_data()
    C_g = calc_coefficients(d, 'Gly', 'n', 'density')
    C_k = calc_coefficients(d, 'MKP', 'n', 'density')

    # if specified difference
    if ratio < 1:
        diff = ratio
        # calculate n
        n = (diff - (C_k[1] - C_g[1])) / (C_k[0] - C_g[0])
    # if specified ratio
    elif ratio >= 1:
        n = - (C_k[1] - ratio * C_g[1]) / (C_k[0] - ratio * C_g[0])

    # calculate densities of fluids
    r_k = C_k[0] * n + C_k[1]
    r_g = C_g[0] * n + C_g[1]

    # volumes of fluid
    v_lock = 0.25 * 0.25 * 0.2
    v_flume = 5 * 0.25 * 0.2

    # calculate wt. % from n
    M_g = calc_coefficients(d, 'Gly', 'n', 'wt.')
    M_k = calc_coefficients(d, 'MKP', 'n', 'wt.')

    # % wt.
    mwt_g = M_g[0] * n + M_g[1]
    mwt_k = M_k[0] * n + M_k[1]

    # absolute mass
    m_g = round(v_flume * r_g * 1000 * mwt_g / 100, 3)
    m_k = round(v_lock * r_k * 1000 * mwt_k / 100, 3)

    # unit costs (gbp/kg)
    ucost_g = 1.1
    ucost_k = 25

    # total costs
    tcost_g = round(ucost_g * m_g, 2)
    tcost_k = round(ucost_k * m_k, 2)
    total = round(tcost_k + tcost_g, 2)

    # solution volume specfic mass (mass of solute per litre solution)
    # mv_g = r_g * mwt_g / 100
    # mv_k = r_k * mwt_k / 100

    # density of water (g/cm^3)
    r_w = 0.9982
    # water volume specific mass (mass of solute per litre water)
    mvw_g = mwt_g / (100 - mwt_g) * r_w
    mvw_k = mwt_k / (100 - mwt_k) * r_w

    print("Density ratio/diff of %s" % ratio)
    print("Requires n = %s" % round(n, 4))
    print("MKP: density = %.4f, total mass = %skg @ %.2f%%mass (%.3fg / 250ml water) " % (r_k, m_k, mwt_k, mvw_k * 1000 / 4))
    print("cost = %sgbp @ %sgbp/kg" % (tcost_k, ucost_k))
    print("Gly: density = %.4f, total mass = %skg @ %.2f%%mass (%.3fg / 250ml water) " % (r_g, m_g, mwt_g, mvw_g * 1000 / 4))
    print("cost = %sgbp @ %sgbp/kg" % (tcost_g, ucost_g))
    print("Total cost is %sgbp" % total)

    return total


def salt(rho, volume=200):
    """Given a density (rho) and (optionally) volume (in litres), calculate
    how much salt is needed.
    """
    d = get_data()
    C_s = calc_coefficients(d, 'NaCl', 'density', 'wt.')
    # given rho, calculate %wt.
    wt = rho * C_s[0] + C_s[1]

    # absolute mass
    m = round(volume * rho * wt / 100, 3)
    # density of water (g/cm^3)
    r_w = 0.9982
    # water volume specific mass (mass of solute per litre water)
    mvw_s = wt / (100 - wt) * r_w
    # solution volume specfic mass (mass of solute per litre solution)
    # mv_s = r_w * wt / 100

    mass_of_scoop = 0.045
    mass_of_level_scoop = 0.425
    level_scoop_salt = mass_of_level_scoop - mass_of_scoop
    no_scoops = round(m / level_scoop_salt, 2)

    new_volume = round((m + r_w * volume) / rho, 0)

    print "You have a volume of {volume}L".format(volume=volume)
    print "To get a density of {rho}, use {mass}kg of salt.".format(rho=rho, mass=m)
    print "The volume after adding salt will be {nv} L".format(nv=new_volume)
    print "That's {rel}kg/L of water".format(rel=round(mvw_s, 3))
    print "Which is about {no_scoops} scoops in total".format(no_scoops=no_scoops)

    return m


def density(n, substance, d=None):
    """Calculate density of a substance at a given value
    of n, using a spline representation.

    The difficulty with this is the extrapolation outside
    the data points given in the crc-data. Setting k=1
    in the spline creation forces our interpolation to be
    linear, giving roughly sane results outside the data
    range.
    """
    if not d:
        d = get_data()
    R = d[substance]['density']
    N = d[substance]['n']
    spline = interpolate.UnivariateSpline(N, R, s=0, k=1)
    density = spline(n)
    return density


def viscosity(n, substance, d=None):
    """Calculate viscosity of a substance at a given value
    of n.

    The variation of viscosity is definitely non-linear.
    We fit a spline to it that goes through all the points.
    """
    if not d:
        d = get_data()
    V = d[substance]['viscosity']
    N = d[substance]['n']
    spline = interpolate.UnivariateSpline(N, V, s=0, k=1)
    visc = spline(n)
    return visc


def S(rc='MKP', r1='NaCl', r2='Gly', n=None, d=None):
    """Calculate the stratification parameter.

    Specify either the densities or the refractive index.

    If only densities specified, rc, r1, r2 are the densities of the
    current, lower and upper layers.

    n is the refractive index. If this is specified, the densities
    of the layers will be worked out for the substances given as the
    density arguments. The corresponding S is returned.
    """
    if not n:
        S = (r1 - r2) / (rc - r2)
    elif n:
        if not d:
            d = get_data()
        rc = density(n, rc, d)
        r1 = density(n, r1, d)
        r2 = density(n, r2, d)
        S = (r1 - r2) / (rc - r2)
    return S


def plot_cost():
    r = np.linspace(1, 1.10)
    c = map(cost, r)
    plt.plot(r, c)
    plt.xlabel(r'ratio of densities, $\rho_c / \rho_0$')
    plt.ylabel(u'Cost (\u00A3)')
    plt.xlim(1, 1.11)
    plt.savefig('cost_vs_density.png')


def plot():
    d = get_data()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # gly_n = d['Gly']['n']
    # nacl_n = d['NaCl']['n']
    # mkp_n = d['MKP']['n']
    # gly_r = d['Gly']['density']
    # nacl_r = d['NaCl']['density']
    # mkp_r = d['MKP']['density']

    # nacl_eta = d['NaCl']['viscosity']
    # sub_gly = d['Gly']['viscosity'][:len(nacl_eta)]
    # deta_gly_nacl = ( sub_gly / nacl_eta )

    # gly_c = np.polyfit(gly_n, gly_r, 1)
    # nacl_c = np.polyfit(nacl_n, nacl_r, 1)
    # mkp_c = np.polyfit(mkp_n, mkp_r, 1)
    gly_c = calc_coefficients(d, 'Gly', 'n', 'density')
    nacl_c = calc_coefficients(d, 'NaCl', 'n', 'density')
    mkp_c = calc_coefficients(d, 'MKP', 'n', 'density')

    n = np.linspace(1.333, 1.35)
    gly_R = n * gly_c[0] + gly_c[1]
    nacl_R = n * nacl_c[0] + nacl_c[1]
    mkp_R = n * mkp_c[0] + mkp_c[1]

    ax.plot(n, mkp_R - gly_R, label=r'$\Delta\rho_{mkp-gly}$')
    ax.plot(n, nacl_R - gly_R, label=r'$\Delta\rho_{nacl-gly}$')
    # ax.plot(gly_n, gly_r, label='Gly')
    # ax.plot(nacl_n, nacl_r, label='NaCl')
    # ax.plot(mkp_n, mkp_r, label='MKP')
    ax.set_xlim(1.333, 1.35)
    ax.set_ylim(0, 0.1)
    ax.set_xlabel('Refractive index, n')
    ax.set_ylabel(r'Density difference, $g/cm^3$')
    ax.grid()
    ax.legend(loc=0)

    # ax2 = ax.twinx()
    # ax2.plot(nacl_n, deta_gly_nacl, label='delta eta')
    # ax2.legend(loc=0)

    plt.show()


def compare_combinations():
    """If we want to set up a two layer system for a gravity current
    to propagate into and we want refractive index matching across
    the entire domain, we are limited in the ratios of density that
    we can choose.

    This function creates a plot comparing various chemical combinations.
    Plotted is the stratification parameter,

    S = (rho_1 - rho_2) / (rho_c - rho_2)

    over the refractive index, n, as this is the control parameter.

    We want the difference between densities of adjacent fluids to
    be at least 0.01 to bring the errors in setting the density in
    check. Lines are only plotted where this is true.

    We are also limited by the solubility of the chemicals in water,
    as this sets the maximum density. Lines are plotted up to this limit.

    If a desirable S and ratio is possible for a given combination,
    we should then check whether the viscosities vary greatly between
    the fluids and seek to minimise this.

    The most important viscosity difference is that between the current
    and the first layer as it is this interface that will be turbulent.

    The viscosity difference between these two layers in percent is
    plotted, with horizontal lines at 5 and 10 %.
    """
    # range of refractive indices to consider
    nlow = 1.3333
    nhi = 1.3600
    N = np.linspace(nlow, nhi)

    combs = [{'rc': 'MKP', 'r1': 'NaCl', 'r2': 'Gly'},
             {'rc': 'MKP', 'r1': 'DKP', 'r2': 'Gly'},
             {'rc': 'MKP', 'r1': 'DKP', 'r2': 'MgCl'},
             # {'rc': 'MKP', 'r1': 'DNP', 'r2': 'Gly'},
             # {'rc': 'MKP', 'r1': 'MNP', 'r2': 'Gly'},
             # {'rc': 'MNP', 'r1': 'DNP', 'r2': 'Gly'},
             {'rc': 'MKP', 'r1': 'MgCl', 'r2': 'Gly'},
             {'rc': 'DKP', 'r1': 'NaCl', 'r2': 'Gly'},
             {'rc': 'DKP', 'r1': 'MgCl', 'r2': 'Gly'}]

    fig = plt.figure()
    axs = fig.add_subplot(1, 1, 1)

    axs.set_ylim(0, 1)
    axs.set_ylabel(r'$S = \frac{\rho_1 - \rho_2}{\rho_c - \rho_2}$')
    axs.set_xlim(nlow, nhi)
    axs.set_xlabel('Refractive index')

    axv = axs.twinx()
    axv.set_xlim(1.3333, 1.3540)
    axv.set_ylim(0, 100)
    axv.set_ylabel('Viscosity difference (%)')

    d = get_data()
    for i, comb in enumerate(combs):
        rc = density(N, comb['rc'], d)
        r1 = density(N, comb['r1'], d)
        r2 = density(N, comb['r2'], d)

        mrc = max_density(comb['rc'], d)
        mr1 = max_density(comb['r1'], d)
        mr2 = max_density(comb['r2'], d)

        R12 = r1 - r2
        Rc1 = rc - r1

        # only where density difference > 0.01
        delta = 0.01
        cond1 = (R12 > delta) & (Rc1 > delta)
        # only where not saturated
        cond2 = (rc < mrc) & (r1 < mr1) & (r2 < mr2)
        cond = cond1 & cond2

        Sn = np.array([S(n=n, d=d, **comb) for n in N])  # can't vectorise keywords!
        Nc = N[np.where(cond)]
        Snc = Sn[np.where(cond)]
        label = "{rc}-{r1}-{r2}".format(rc=comb['rc'], r1=comb['r1'], r2=comb['r2'])
        axs.plot(N, Sn, 'k', alpha=0.1)
        axs.plot(Nc, Snc, label=label)

        V = np.abs(1 - viscosity(N, comb['r1'], d) / viscosity(N, comb['rc'], d)) * 100
        Vc = V[np.where(cond)]

        axv.plot(Nc, Vc)

    leg = axs.legend(loc='upper left', ncol=1)
    leg.legendPatch.set_alpha(0.3)  # legend transparency
    plt.setp(leg.get_texts(), fontsize='small')
    axv.axhline(5, color='k', linestyle='--')
    axv.axhline(10, color='k', linestyle='--')

    fig.savefig('Svn-visc.png')
    return fig


def compare_substances(n=1.3450, dn=0, step=5):
    """Plot substances in (density, viscosity) for a given refractive
    index.

    TODO:
        Optionally plot for two refractive indices to get an idea of
        behaviour with n for each substance. Connect up two points with
        a line.
    """
    d = get_data()
    substances = d.keys()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    N = np.linspace(n, n + dn, step)

    for sub in substances:
        R = density(N, sub, d)
        V = viscosity(N, sub, d)
        ax.plot(R, V, 'o')
        ax.annotate(sub, xy=(R[0], V[0]),
                    xytext=(-20, 5), textcoords='offset points')

    ax.set_xlim(1.00, 1.10)
    ax.set_xlabel('Density')
    ax.set_ylim(0.95, 1.50)
    ax.set_ylabel('Viscosity')
    title = 'Comparison of aqueous solutions at n={n}'.format(n=n)
    ax.set_title(title)

    fig.savefig('chemical-comparison.png')

    return fig


if __name__ == '__main__':
    cost(float(argv[1]))
