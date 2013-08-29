Refractive index matching calculator
------------------------------------

### Aim ###

Perform an experiment with two fluid phases of different densities
that have matched refractive indices.

### Why? ###

Accurate Particle Imaging Velocimetry (PIV) requires a homogeneous
refractive index in the field of view. Refractive index (RI) varies
with density and is a function of the substance used to vary the
density.

In a two density system it is necessary to match the RI of the two
fluids to have a homogeneous RI field. This means using two different
substances with difference RI / density relations in order to match
the RI whilst maintaining density separation.

For example, one might use Glycerol and Monopotassium Phosphate as a
pair of substances. Saline and Ethyl alcohol also work. In fact
there are lots of valid combinations of substances (see the optical
properties of aqueous solutions section in the CRC handbook), but we
are constrained by cost and availability.

### Usage ###

Initialise a run with a given target density ratio and the
substances to be used, and (optionally) the volumes of the two fluid
phases:

    ```python
    run = RIMatched(density=1.05, V1=100, V2=200, sub1='MKP', sub2='Gly')
    ```

You can refer to substance 1 as any of `run.sub1`, `run.mkp`, or any
name that you supply at instantiation with the named argument
`name1`.

Output the quantity of each substance needed to achive the target
density ratio:

    ```python
    run.quantities
    ```

How much is this all going to cost?

    ```python
    run.total_cost
    ```

Given the density that you have mixed one phase up to, how much more
substance (or water, if you've overshot the density) is needed to
make the target?

    ```python
    run.substance.how_much_more(density_measurement=what_you_measured)
    ```

Equally, given a refractive index measurement, what do you need to
add to get to the target?

    ```python
    run.substance.how_much_more(ri_measurement=what_you_measured)
    ```


You're nearly ready to go, but there is a temperature difference
between the two phases. Given a RI measurement of one phase, what
measured RI do you need to shoot for for the other phase, and how do
you get there?

TODO: spec temp calibration and implement
TODO: use full dn/dt data
TODO: codify sequence of operations for typical lab expt.

You've matched the RI and want to log the final measurements of
density, RI and temperature:

TODO: implement logging
    ```python
    run.log(density1=1.06, density2=1.04, temperature1=.. 
    ```

You want to append the log entry to a log file with the rest of your
experimental logs:

TODO: implement logging
    ```python
    run.save(logfile)
    ```
