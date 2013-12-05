from __future__ import division

import argparse

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import interpolate

# density of water (g/cm^3)
density_water = 0.9982

# unit costs (gbp/kg)
unit_cost = {'MKP': 2.8,
             'Gly': 1.1}

# tank volume in litres
v_lock = 0.25 * 0.25 * 0.2 * 1000
v_flume = 5 * 0.25 * 0.2 * 1000

# masses in kg
mass_of_scoop = 0.045
mass_of_level_scoop = 0.425
level_scoop_salt = mass_of_level_scoop - mass_of_scoop

# approx dn/dt for glycerol and mkp
# TODO: update these with more data
# TODO: separate data and calc and move calc inside AqueousSolution
dndt = {'Gly': (1.3370 - 1.3381) / (19.0 - 8.8),
        'MKP': (1.3368 - 1.3380) / (19.8 - 8.8)}


class AqueousSolution(object):
    """A mixture of a given chemical that is used as a solute and
    water as a solvent, mixed to a given refractive index or density
    with a given initial volume of water.

    Various methods for calculating properties of the solute and the
    solution.
    """
    def __init__(self, ref, n=None, density=None, volume=1, temperature=20):
        """
        Inputs: ref     - string, e.g. 'MKP', 'Gly'
                n       - float, refractive index
                density - float, density (g/(cm^3))
                volume  - float, volume of water in litres
                temperature - float, temperature of bulk solution in C
                              (used for refractive index)

        n and density default to a zero concentration solution. If
        you set either one of them, the other will be (re)calculated.

        You cannot set both n and density for a specific substance.
        """
        self.ref = ref
        self.data = get_data()[self.ref]
        self.volume = volume
        self.temperature = temperature

        self.dndt = dndt[self.ref]

        if not n and not density:
            self.target_density = density_water
        elif not n:
            self.target_density = density
        elif not density:
            self.target_n = n
        else:
            raise UserWarning("You can't constrain both "
                              "refractive index and density!")

    @property
    def target_density(self):
        return self._density

    @target_density.setter
    def target_density(self, rho):
        """Set the density, raising an error if it exceeds the
        maximum possible density for the solution."""
        md = self.max_density
        if rho > md:
            msg = ('Not possible to set density to {rho:.4f}. Maximum density '
                   'for {sub} is {md:.4f}').format(rho=rho,
                                                   sub=self.ref,
                                                   md=md)
            raise UserWarning(msg)
        self._density = rho
        # calculate the refractive index we would have at T=20
        # (at which CRC data is measured)
        n_TC = self.n(rho).flatten()[0]
        # convert this to the r.i. we would have at solution
        # temperature
        self._n = n_TC + self.dndt * (self.temperature - 20)

    @property
    def target_n(self):
        """The refractive index that the bulk solution should
        be at at 20C, correcting for temperature (i.e. if T!=20C)
        """
        return self._n

    @property
    def target_n_tc(self):
        """The refractive index that the bulk solution should
        be at if T=20C.
        """
        return self.calc_n_sample(T_sample=20)

    @target_n.setter
    def target_n(self, n):
        self._n = n
        # calculate the refractive index that we would have at T=20C
        n_TC = self.calc_n_sample(T_sample=20)
        # flatten()[0] is to extract value from a 0d array
        self._density = self.density(n_TC).flatten()[0]

    def calc_n_sample(self, T_sample=20):
        """The refractive index that will be measured at a given
        temperature, T_sample.
        """
        return (self.temperature - T_sample) * self.dndt + self.target_n

    @property
    def solution_volume(self):
        """Volume of solution after mixing in solute to initial
        volume of water.
        """
        rho_soln = self.target_density
        rho_water = density_water
        v_water = self.volume
        m_solu = self.absolute_mass

        v_soln = (m_solu + rho_water * v_water) / rho_soln
        return v_soln

    def set_volume(self, V):
        """Set the volume of the substance."""
        self.volume = V

    def calc_coefficients(self, x, y):
        """Assuming a straight line fit, calculate (m,c) for specific
        values of x, y (density, n, etc.) for a specific substance.

        x         - string, e.g. 'n'
        y         - string, e.g. 'density'

        i.e. fit y = m * x + c and return (m, c)
        """
        m, c = np.polyfit(self.data[x], self.data[y], 1)
        return m, c

    def density(self, n):
        """Calculate density of a substance at a given value
        of n, using a spline representation.

        The difficulty with this is the extrapolation outside
        the data points given in the crc-data. Setting k=1
        in the spline creation forces our interpolation to be
        linear, giving roughly sane results outside the data
        range.
        """
        d = self.data
        R = d['density']
        N = d['n']
        spline = interpolate.UnivariateSpline(N, R, s=0, k=1)
        density = spline(n)
        return density

    def n(self, density):
        """Calculate refractive index of a substance at a given
        value of density, using a spline representation.
        """
        d = self.data
        R = d['density']
        N = d['n']
        spline = interpolate.UnivariateSpline(R, N, s=0, k=1)
        n = spline(density)
        return n

    @property
    def max_density(self):
        """Calculate the maximum density that a substance can be
        mixed to in water.
        """
        m, c = self.calc_coefficients('wt.', 'density')
        d = self.data
        sol = d['solubility']
        wt = (sol / (sol + 1000)) * 100
        max_density = m * wt + c
        return max_density

    def viscosity(self, n):
        """Calculate viscosity of a substance at a given value
        of n.

        The variation of viscosity is definitely non-linear.
        We fit a spline to it that goes through all the points.
        """
        d = self.data
        V = d['viscosity']
        N = d['n']
        spline = interpolate.UnivariateSpline(N, V, s=0, k=1)
        visc = spline(n)
        return visc

    @property
    def target_percent_weight(self):
        """Calculate the % weight of solution of given substance, to
        achieve the target refractive index.

        % weight is (mass of solute) / (total mass of solution) * 100.
        """
        M = self.calc_coefficients('n', 'wt.')
        mwt = M[0] * self.target_n + M[1]
        return mwt

    @property
    def absolute_mass(self):
        """Calculate the absolute mass of given substance needed to
        achieve the target density.
        """
        m_sub = self.specific_mass * self.volume
        return m_sub

    @property
    def specific_mass(self):
        """Calculate specific mass of given substance (kg / litre of water)."""
        # water volume specific mass (mass of solute per litre water)
        mwt_sub = self.target_percent_weight
        mvw_sub = mwt_sub / (100 - mwt_sub) * density_water
        return mvw_sub

    @property
    def specific_mass_g(self):
        """Specific mass in g / L."""
        return self.specific_mass * 1000

    @property
    def specific_mass_solution(self):
        """Calculate specific mass of given substance (kg / litre of
        solution).
        """
        # solution volume specfic mass (mass of solute per litre solution)
        mwt_sub = self.target_percent_weight
        r_sub = self.target_density
        mv_sub = r_sub * mwt_sub / 100
        return mv_sub

    @property
    def specific_mass_solution_g(self):
        """Specific mass in g / L."""
        return self.specific_mass_solution * 1000

    @property
    def unit_cost(self):
        """Return the cost of a substance in GBP / kg."""
        unit_cost_sub = unit_cost.get(self.ref, 0)
        return unit_cost_sub

    @property
    def cost(self):
        """Return the total cost of a particular substance, in GBP."""
        unit_cost_sub = self.unit_cost
        mass_sub = self.absolute_mass
        tcost_sub = round(unit_cost_sub * mass_sub, 2)
        return tcost_sub

    def n_sensitivity(self, n, volume=None, dn=0.0001):
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
        volume = volume or self.volume
        rho = self.density(n)
        # total mass
        M = volume * rho * 1000
        # variation of n with % wt.
        m, c = self.calc_coefficients('wt.', 'n')
        dwt = dn / m
        dm = M * dwt / 100

        return dm

    @property
    def no_scoops(self):
        """How many lab scoops of solute is this?"""
        m = self.absolute_mass
        no_scoops = round(m / level_scoop_salt, 2)
        return no_scoops

    @property
    def instructions(self):
        """Return a string of instructions for mixing up given substance."""
        ins_str = u"""
        {s.ref}: density = {s.target_density:.4f} g / (cm)^3,
             volume = {s.volume} L,
             temperature = {s.temperature} C,
             total mass = {s.absolute_mass} kg
                          @ {s.target_percent_weight:.2f} %mass
             ({s.specific_mass_g:.3f}g / L water)
             solution volume = {s.solution_volume:.4f} L
             # level scoops = {s.no_scoops}
             cost = {s.cost} gbp @ {s.unit_cost} gbp/kg
             refractive index at T=20: {s.target_n_tc}
        """.format(s=self)
        return ins_str

    def how_much_more(self, density_measurement=None, ri_measurement=None,
                      density_tolerance=None, ri_tolerance=None):
        """Given a density measurement (g / (cm)^3) or a refractive
        index measurement, how much more solute or water is needed
        to make the target density?
        """
        if not density_measurement and not ri_measurement:
            raise UserWarning("Need to specify either density or RI")
        # measured density
        rho_m = density_measurement or self.density(ri_measurement)
        # target density
        rho_t = self.target_density
        # water density
        rho_0 = density_water

        ri_difference = ri_measurement - self.target_n
        density_difference = rho_m - rho_t

        if (np.abs(density_difference) < density_tolerance) \
           or (np.abs(ri_difference) < ri_tolerance):
            # within measurement error, add nothing
            # < None is always False
            return 'Matched, add nothing'

        elif rho_m < rho_t:
            # add solute
            sub_m = AqueousSolution(self.ref,
                                    density=rho_m,
                                    volume=self.volume)
            how_much_solute = self.absolute_mass - sub_m.absolute_mass
            return (self.ref, how_much_solute, 'kg',
                    'dn=', round(ri_difference, 5))

        elif rho_m > rho_t:
            # add water
            how_much_water = self.solution_volume * (rho_m - rho_t) \
                / (rho_t - rho_0)
            return ('Water', how_much_water, 'L',
                    'dn=', round(ri_difference, 5))


class RIMatched(object):
    """Represents a two phase fluid experiement in which two
    substances are used to mix the fluids to a given density ratio
    whilst maintaining equal refractive indices.
    """
    def __init__(self, density_ratio, sub1='MKP', sub2='Gly',
                 v1=v_lock, v2=v_flume, t1=20, t2=20,
                 name1='lock', name2='tank', density_floor=1.00):
        """
        Inputs: density_ratio - float, the target density ratio,
                                i.e. rho_1 / rho_2
                sub1    - string, substance 1, default 'MKP'(Monopotassium
                          Phosphate)
                sub2    - string, substance 2, default 'Gly' (Glycerol)
                v1      - float, Volume of substance 1 in experiment (litres)
                v2      - float, Volume of substance 2 in experiment (litres)
                t1      - float, operating temperature of substance 1 (C)
                t2      - float, operating temperature of substance 2 (C)

                density_floor - float, default 1.0, minimum density of either
                                of the layers
        """
        self.ratio = density_ratio
        self.density_floor = density_floor
        # Each substance is a AqueousSolution
        self.sub1 = AqueousSolution(sub1, volume=v1, temperature=t1)
        self.sub2 = AqueousSolution(sub2, volume=v2, temperature=t2)
        # different aliases for substances
        setattr(self, sub1.lower(), self.sub1)
        setattr(self, sub2.lower(), self.sub2)
        self.name1 = name1
        self.name2 = name2
        setattr(self, self.name1, self.sub1)
        setattr(self, self.name2, self.sub2)

        # calculate the required refractive index and set it on the
        # substances
        n = self.n_matched
        self.sub1.target_n = n
        self.sub2.target_n = n

        self.V1 = v1
        self.V2 = v2

        # volumes in m^3
        self.V1m3 = v1 / 1000
        self.V2m3 = v2 / 1000

    @property
    def n_matched(self):
        """Calculate the refractive index needed to achieve the
        target density ratio.
        """
        # functions that calculate density as a function of n
        f1 = self.sub1.density
        f2 = self.sub2.density

        # find intersection. monotonic, so look for f1 - f2 = 0
        # FIXME: for density floors > 1, you start to limit the
        # range of n_matched densities. account for this.
        # Basically, the function is not smooth - not sure if
        # this is actually a problem.
        def f(n):
            d1 = f1(n)
            d2 = f2(n)
            if d2 < self.density_floor:
                d2 = self.density_floor
            return d1 / d2 - self.ratio

        n = scipy.optimize.bisect(f, 1.3, 1.5)

        return n

    @property
    def total_cost(self):
        """Return total cost of run in GBP."""
        tcost_g = self.sub2.cost
        tcost_k = self.sub1.cost

        total = round(tcost_k + tcost_g, 2)
        return total

    def total_cost_instructions(self):
        """Takes desired density difference or ratio as input and
        calculates the cost of a run based on unit costs of glycerol
        and mkp.

        Outputs: string, instructions including cost of run,
                 quantities of gly/mkp, required n.
        """
        n = self.n_matched

        ins = []
        ins.append("Density ratio of {ratio}".format(ratio=self.ratio))
        ins.append("Requires n = {n}".format(n=round(n, 4)))
        ins.append(self.sub1.instructions)
        ins.append(self.sub2.instructions)
        ins.append("Total cost is {total}gbp".format(total=self.total_cost))

        instructions = '\n'.join(ins)

        return instructions

    @property
    def quantities(self):
        print(self.total_cost_instructions())


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
    solubilities = [('Gly', 9999),      # units g substance / L water
                    ('NaCl', 359),      # (wikipedia)
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


def S(rc='MKP', r1='NaCl', r2='Gly', n=None):
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
        rc = AqueousSolution(rc).density(n)
        r1 = AqueousSolution(r1).density(n)
        r2 = AqueousSolution(r2).density(n)
        S = (r1 - r2) / (rc - r2)
    return S


def plot_cost():
    r = np.linspace(1, 1.10)
    Exps = map(RIMatched, r)
    c = [e.total_cost for e in Exps]
    plt.plot(r, c)
    plt.xlabel(r'ratio of densities, $\rho_c / \rho_0$')
    plt.ylabel(u'Cost (\u00A3)')
    plt.xlim(1, 1.11)
    plt.savefig('cost_vs_density.png')


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # nacl_eta = d['NaCl']['viscosity']
    # sub_gly = d['Gly']['viscosity'][:len(nacl_eta)]
    # deta_gly_nacl = ( sub_gly / nacl_eta )

    gly = AqueousSolution('Gly')
    nacl = AqueousSolution('NaCl')
    mkp = AqueousSolution('MKP')

    n = np.linspace(1.333, 1.35)
    gly_R = gly.density(n)
    nacl_R = nacl.density(n)
    mkp_R = mkp.density(n)

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

    for i, comb in enumerate(combs):
        sc = AqueousSolution(comb['rc'])
        s1 = AqueousSolution(comb['r1'])
        s2 = AqueousSolution(comb['r2'])

        rc = sc.density(N)
        r1 = s1.density(N)
        r2 = s2.density(N)

        mrc = sc.max_density
        mr1 = s1.max_density
        mr2 = s2.max_density

        R12 = r1 - r2
        Rc1 = rc - r1

        # only where density difference > 0.01
        delta = 0.01
        cond1 = (R12 > delta) & (Rc1 > delta)
        # only where not saturated
        cond2 = (rc < mrc) & (r1 < mr1) & (r2 < mr2)
        cond = cond1 & cond2

        Sn = np.array([S(n=n, **comb) for n in N])  # can't vectorise keywords!
        Nc = N[np.where(cond)]
        Snc = Sn[np.where(cond)]
        label = "{sc.ref}-{s1.ref}-{s2.ref}".format(sc=sc, s1=s1, s2=s2)
        axs.plot(N, Sn, 'k', alpha=0.1)
        axs.plot(Nc, Snc, label=label)

        V = np.abs(1 - s1.viscosity(N) / sc.viscosity(N)) * 100
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
        s = AqueousSolution(sub)
        R = s.density(N)
        V = s.viscosity(N)
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
    parser = argparse.ArgumentParser()
    parser.add_argument('density',
                        help=("Density ratio to calculate chemical"
                              "quantities for."),
                        type=float,
                        nargs='?',
                        default=None)
    parser.add_argument('--chem',
                        help=("Specific substances to calculate quantities "
                              "for to achieve density given."),
                        nargs='*')
    parser.add_argument('--list_chems',
                        help=("List possible chemical names"),
                        action='store_true')
    parser.add_argument('--volume',
                        help=("Volume of solution to calculate for"
                              "to achieve density given."),
                        type=float,
                        default=200)
    parser.add_argument('--ri',
                        help=("Refractive index of solution"),
                        nargs=1,
                        type=float)
    parser.add_argument('--v1',
                        help=("Volume of substance 1"),
                        nargs='?',
                        default=v_lock,
                        type=float)
    parser.add_argument('--v2',
                        help=("Volume of substance 2"),
                        nargs='?',
                        default=v_flume,
                        type=float)
    parser.add_argument('--t1',
                        help=("Temperature of substance 1"),
                        nargs='?',
                        default=20.,
                        type=float)
    parser.add_argument('--t2',
                        help=("Temperature of substance 2"),
                        nargs='?',
                        default=20.,
                        type=float)
    args = parser.parse_args()

    if args.list_chems:
        d = get_data()
        print("\nPossible substance names:\n")
        print(d.keys())
        print("")

    elif args.chem:
        for sub in args.chem:
            s = AqueousSolution(sub, density=args.density,
                                volume=args.volume,
                                n=args.ri)
            print(s.instructions)

    else:
        r = RIMatched(density_ratio=args.density,
                      v1=args.v1, v2=args.v2,
                      t1=args.t1, t2=args.t2)
        print(r.total_cost_instructions())
