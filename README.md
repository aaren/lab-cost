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
    run = RIMatched(density=1.05, V1=100, V2=200, substance1='MKP', substance2='GLY')
    ```

Output the quantity of each substance needed to achive the target
density ratio:

    ```python
    run.quantities
    ```

Given the density that you have mixed one phase up to, how much more
substance (or water, if you've overshot the density) is needed to
make the target?

    ```python
    run.how_much_more(density_measurement, substance)
    ```

You're nearly ready to go, but there is a temperature difference
between the two phases. Given a RI measurement of one phase, what
measured RI do you need to shoot for for the other phase, and how do
you get there?



You've matched the RI and want to log the final measurements of
density, RI and temperature:

    ```python
    run.log(density1=1.06, density2=1.04, temperature1=.. 
    ```

You want to append the log entry to a log file with the rest of your
experimental logs:

    ```python
    run.save(logfile)
    ```
