# LFQPO comptonization 

A numerical code for computing Compton scattered X-ray time-dependent spectra arising from a hot inner flow which undergoes Lense-Thirring precession. 

### What does it do
This Monte-Carlo code is supposed to compute the Compton scattered X-ray flux arising from a hot inner flow which undergoes Lense-Thirring precession. The hot flow intercepts seed photons from an outer truncated thin disk. A fraction of the Comptonized photons illuminates back the disk and the reflected/reprocessed photons contribute to the observed spectrum. The total spectrum, including disk thermal emission, hot flow Comptonization, and disk reflection, is modelled within the framework of general relativity, taking light-bending and gravitational/Doppler shift into account. The simulations are performed in the context of the Lense-Thirring precession model for the low-frequency quasi-periodic oscillations observed in X-ray flux of accreting black holes, so the inner flow is assumed to precess leading to a periodic modulation of the emitted radiation.

More information on this code and the simulation results can be found in the scientific paper [arXiv:1801.04028](https://arxiv.org/abs/1801.04028) ([ADS link](http://adsabs.harvard.edu/abs/2018arXiv180104028Y))

### Usage

##### Prerequisities
The code requires MPI and GSL library.

##### Code parts
* absorp.dat (database file of opacity)
* mtc_incl_def.c (definitions)
* mtc_incl_wedge.c (wedge geometry)
* quadrat.c (ODE integration)
* mtc_incl_code.c (code body)
* sim5lib.c ([SIM5 library](https://github.com/mbursa/sim5))
* te.dat (corona temperature according to [arXiv:1801.04028](https://arxiv.org/abs/1801.04028))

##### Step 1 
Compile and run `change_inp.c` to get a series of "mtc__inp.0000XX", where XX  
spans from 00 to 63, corresponding to the 64 phase bins over one complete precession cycle. 
```
$ gcc change_inp.c -o change_inp 
$ ./change_inp
```
The main input parameters of the code, e.g., black hole spin, truncation radius   
and precession angles, are hard-coded in `change_inp.c`, so you need to edit the source to change them.

##### Step 2
Compile and run `sphere_table.c` to get "sphere_table.dat" file (the size of that file will be ~30MB).
```
$ gcc sphere_table.c -lm -lgsl -lgslcblas -o sphere_table
$ ./sphere_table
```
##### Step 3
Compile and run `photon_table.c` to get a series of "table0000XX.dat", where XX spans from 0 to 63, corresponding to the 64 phase bins over one complete precession cycle. 
```
$ gcc photon_table.c -lm -lgsl -lgslcblas -g -o photon_table
$ for i in `seq -w 0 63`; do ./photon_table $i; done
```
Note that the size of each "table0000XX.dat" file will be ~3GB and moreover that this step requires ~3GB of free memory space (RAM) of the computer. In most cases, it will fail on a normal laptop.

##### Step 4
Compile and run `lfqpo_mpi.c` to get the final results (you need OpenMPI).
```
$ mpicc lfqpo_mpi.c -lm -lgsl -lgslcblas -g -o lfqpo
$ mpirun -np 4 ./lfqpo
```

The products of the calculation will consist of:
* accretion disk spectrum: inpsp0000*.dat
* comptonization spectrum: mcomp0000*.dat
* disk reflection spectrum: refls0000*.dat
* Fe Ka emission line: lines0000*.dat

The structure of the spectrum file is:
* E = photon energy
* N = photon counts
* F_E  = normalized spectrum flux

| First 4 columns are for 4 azimuths of observer | then repeat for the next inclination |
| ------ | ------ |
| (E,N1,F_E1), (E,N2,F_E2), (E,N3,F_E3) , (E,N4,F_E4) | ...................................................... |

If you want to plot the total spectrum, please add F_E for disk, Comptonization, reflection and line, rather than adding N.

Should you have any question, please contact Dr. Bei You (youbeiyb@gmail.com).