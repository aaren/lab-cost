from __future__ import division

from sys import argv

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# density of water (g/cm^3)
density_water = 0.9982

# unit costs (gbp/kg)
unit_cost = {'MKP': 2.8,
             'Gly': 1.1}

# tank volume in litres
v_lock = 0.25 * 0.25 * 0.2 * 1000
v_flume = 5 * 0.25 * 0.2 * 1000


class RIMatched(object):
    """Represents a two phase fluid experiement in which two
    substances are used to mix the fluids to a given density ratio
    whilst maintaining equal refractive indices.
    """
    def __init__(self, density_ratio, sub1='MKP', sub2='Gly',
                 V1=v_lock, V2=v_flume):
        """
        Inputs: density_ratio - float, the target density ratio,
                                i.e. rho_1 / rho_2
                sub1    - string, substance 1, default 'MKP'(Monopotassium
                          Phosphate)
                sub2    - string, substance 2, default 'Gly' (Glycerol)
                V1      - float, Volume of substance 1 in experiment (litres)
                V2      - float, Volume of substance 2 in experiment (litres)
        """
        self.ratio = density_ratio
        self.sub1 = sub1
        self.sub2 = sub2
        self.V1 = V1
        self.V2 = V2

        # volumes in m^3
        self.V1m3 = V1 / 1000
        self.V2m3 = V2 / 1000

        self.data = get_data()

    def volume(self, substance):
        """Return the volume (in litres) of the given substance.

        1L = 0.001 m^3 so this result needs to be divided by 1000 to
        get SI units.
        """
        if substance == self.sub1:
            return self.V1
        elif substance == self.sub2:
            return self.V2
        else:
            self.no_substance_warning()

    @property
    def n_matched(self):
        """Calculate the refractive index needed to achive the
        target density ratio.
        """
        C_2 = calc_coefficients(self.data, self.sub2, 'n', 'density')
        C_1 = calc_coefficients(self.data, self.sub1, 'n', 'density')
        n = - (C_1[1] - self.ratio * C_2[1]) / (C_1[0] - self.ratio * C_2[0])
        return n

    def target_density(self, substance):
        """Calculate the density of given substance needed to achive the
        target density ratio.

        Returns the target density in units of g / (cm)^3.

        1 g / (cm)^3 = 1000 kg / m^3, so this value needs to be
        multiplied by 1000 for it to be in SI units.
        """
        C = calc_coefficients(self.data, substance, 'n', 'density')
        n = self.n_matched
        r = C[0] * n + C[1]
        return r

    @property
    def density1(self):
        """Target density for substance 1."""
        return self.target_density(self.sub1)

    @property
    def density2(self):
        """Target density for substance 2."""
        return self.target_density(self.sub2)

    def target_percent_weight(self, substance):
        """Calculate the % weight of solution of given substance, to
        achieve the target density.

        % weight is (mass of solute) / (total mass of solution) * 100.
        """
        M = calc_coefficients(self.data, substance, 'n', 'wt.')
        mwt = M[0] * self.n_matched + M[1]
        return mwt

    @property
    def percent_weight1(self):
        """Target percent weight for substance 1."""
        return self.target_percent_weight(self.sub1)

    @property
    def percent_weight2(self):
        """Target percent weight for substance 2."""
        return self.target_percent_weight(self.sub2)

    def absolute_mass(self, substance):
        """Calculate the absolute mass of given substance needed to
        achieve the density ratio.

        Returns the mass in kg to 3 decimal places.
        """
        # density of substance in SI units
        r_sub = self.target_density(substance) * 1000
        # percent weight of substance
        mwt_sub = self.target_percent_weight(substance)
        # volume of substance in SI units
        v_sub = self.volume(substance) / 1000

        m_sub = round(v_sub * r_sub * mwt_sub / 100, 3)
        return m_sub

    def absolute_mass1(self):
        """Absolute mass for substance 1."""
        return self.absolute_mass(self.sub1)

    def absolute_mass2(self):
        """Absolute mass for substance 2."""
        return self.absolute_mass(self.sub2)

    def specific_mass(self, substance):
        """Calculate specific mass of given substance (per litre of water)."""
        # water volume specific mass (mass of solute per litre water)
        mwt_sub = self.target_percent_weight(substance)
        mvw_sub = mwt_sub / (100 - mwt_sub) * density_water
        return mvw_sub

    def specific_mass_solution(self, substance):
        """Calculate specific mass of given substance (per litre of solution)."""
        # solution volume specfic mass (mass of solute per litre solution)
        mwt_sub = self.target_percent_weight(substance)
        r_sub = self.target_density(substance)
        mv_sub = r_sub * mwt_sub / 100
        return mv_sub

    def unit_cost(self, substance):
        """Return the cost of a substance in GBP / kg."""
        unit_cost_sub = unit_cost.get(substance, None)
        if not unit_cost_sub:
            self.no_substance_warning()
        return unit_cost_sub

    @staticmethod
    def no_substance_warning(substance):
        msg = ('{substance} not one specified'
                'for this run').format(substance=substance)
        raise UserWarning(msg)

    def cost(self, substance):
        """Return the total cost of a particular substance, in GBP."""
        unit_cost_sub = self.unit_cost(substance)
        mass_sub = self.absolute_mass(substance)
        tcost_sub = round(unit_cost_sub * mass_sub, 2)
        return tcost_sub

    @property
    def total_cost(self):
        """Return total cost of run in GBP."""
        tcost_g = self.cost(self.sub2)
        tcost_k = self.cost(self.sub1)

        total = round(tcost_k + tcost_g, 2)
        return total

    def instructions(self, substance):
        """Return a string of instructions for mixing up given substance."""
        # TODO: consider using attributes of substance in string
        # e.g. sub=substance.name, rho=substance.density, ...
        r_sub = self.target_density(substance)
        v_sub = self.volume(substance)
        mass_sub = self.absolute_mass(substance)
        mwt_sub = self.target_percent_weight(substance)
        mvw_sub = self.specific_mass(substance)
        cost_sub = self.cost(substance)
        ucost_sub = self.unit_cost(substance)
        ins_str = """
        {sub}: density = {rho:.4f} g / (cm)^3,
             volume = {volume} L,
             total mass = {mass}kg @ {mwt:.2f}%mass
             ({mvw:.3f}g / L water)
             cost = {cost}gbp @ {ucost}gbp/kg
        """.format(sub=substance, volume=v_sub, rho=r_sub, mass=mass_sub, mwt=mwt_sub,
                    mvw=mvw_sub * 1000, cost=cost_sub, ucost=ucost_sub)
        return ins_str

    def total_cost_instructions(self):
        """Takes desired density difference or ratio as input and
        calculates the cost of a run based on unit costs of glycerol
        and mkp.

        Outputs: string, instructions including cost of run,
                 quantities of gly/mkp, required n.
        """
        ratio = self.ratio
        d = self.data
        n = self.n_matched

        # calculate densities of fluids
        r_k = self.density1
        r_g = self.density2

        # volumes of fluid in cubic metres
        v_lock = self.V1m3
        v_flume = self.V2m3

        # % wt.
        mwt_g = self.percent_weight2
        mwt_k = self.percent_weight1

        # absolute mass
        m_g = self.absolute_mass2
        m_k = self.absolute_mass1

        # total costs
        tcost_g = self.cost(self.sub2)
        tcost_k = self.cost(self.sub1)

        total = self.total_cost

        # water volume specific mass (mass of solute per litre water)
        mvw_g = self.specific_mass(self.sub2)
        mvw_k = self.specific_mass(self.sub1)

        ins = []
        ins.append("Density ratio of {ratio}".format(ratio=ratio))
        ins.append("Requires n = {n}".format(n=round(n, 4)))
        ins.append(self.instructions(self.sub1))
        ins.append(self.instructions(self.sub2))
        ins.append("Total cost is {total}gbp".format(total=total))

        instructions = '\n'.join(ins)

        return instructions


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
                  ('KBr', 159, -1)]
    solubilities = [('Gly', 9999),      # units g substance / L water. (wikipedia)
                    ('NaCl', 359),
                    ('MKP', 220),
                    ('DKP', 150),
                    ('KCl', 281),
                    ('MNP', 599),
                    ('DNP', 77),
                    ('MgCl', 543),
                    ('LiCl', 830),
                    ('KBr', 620)]
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

    substance - string, e.g. 'MKP'
    x         - string, e.g. 'n'
    y         - string, e.g. 'density'

    i.e. fit y = m * x + c and return (m, c)
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


# very approximately, the dn / dt for glycerol is
dndt = {}
dndt['gly'] = (1.3370 - 1.3381) / (19.0 - 8.8)
# approx dndt for mkp
dndt['mkp'] = (1.3368 - 1.3380) / (19.8 - 8.8)


def ri(chem, t_real, n_sample, t_sample):
    """Calculate the real refractive index of a fluid, given the
    sample refractive index and temperature (n_sample, t_sample)
    and the temperature of the fluid in situ (t_real).
    """
    # the ar200 measures the ri at some T, which will not be the
    # same as the initial sample T.
    n_real = (t_real - t_sample) * dndt[chem] + n_sample
    return n_real


def set_ri_lock(t_lock, t_lock_sample, n_tank):
    # we want to set the lock ri to be the same as the tank ri,
    # n_tank
    # we can only measure n_lock_sample, so find out what this
    # should be for n_lock to equal n_tank
    n_lock = n_tank
    n_lock_sample = (t_lock_sample - t_lock) * dndt['mkp'] + n_lock
    return n_lock_sample


if __name__ == '__main__':
    ratio = float(argv[1])
    r = RIMatched(density_ratio=ratio)
    print(r.total_cost_instructions())
