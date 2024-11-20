Python codes for analyzing atom dump files from LAMMPS runs
All written to run under Python3
Each has a header detailing the program's function, input, and output. 
Each also has parameters that must be adjusted to match parameters used in the simulation run.
Written for me, not you - so examine code to get all details.

calcdisp8.py

> Takes as input the dump and intemp files created when an existing film is thermocycled past its thermal relaxation temperature twice.
> Outputs file of atom mobility and density as a function of temperature for first and second thermocycle.
> Code has to identify deposition, start, and stop temperature

calcdispTT6.py

>Same output as previous but for film held at constant temperature above Tg
>meant to look at dynamics of surface melting into bulk

Cluster_sim.py

>Calculate cluster size distribution for randomly placed 'atoms' in 18x18x35 lattice
> outputs cluster distribution as a function of 'atom' concentration

contactlist.py

>Calculates the average number of contacts for each atom size
>Using dump file containing atoms with no kinetic energy 
>Assumes Lennard-Jones potential between atoms varying by atom type 
 (have to match to values used in MDS)
>Contact if separation less than 2^(⅙)*sigma
>Code specifies z-range for which statistics calculated
>Outputs contact numbers for every combination of atom types
 
expnsncoeff3.py

>Calculatesfilm density at a set of temperatures using files created by LAMMPS cooldown simulation (in.cool-pe-dens)

Filmdens.py

>calculates film density vs distance from growing surface averaged over entire film growth

hyper_uniform_1.py

>Calculates thedistribution of relative sphere sizes of a grown film
>both the number at a distance and the number within a distance
>normalized by the number of atoms studied

TopDensGrad2.py

>Analyze density and size distribution for swapping in >>previously grown<<< film
>Gives both average density and type distribution vs distance from growing surface

atom-atom-sep.py
>Calculate the distribution of atom-atom separations relative to defined diameter
>meant for zero temperature films
