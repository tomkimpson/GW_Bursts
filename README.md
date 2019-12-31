
# Gravitational Burst Radiation from pulsars around massive BHs

This code calculates the waveforms and SNRs of gravitational burst radiation from a pulsar orbiting a massive BH. We use a [numerical kludge method](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.75.024005) where the relativistic spin-orbital dynamics are specified by the Mathsisson Papetrou Dixon formulation (see the [SpinCurvature git repo](https://github.com/tomkimpson/SpinCurvatureDynamics)). 

The orbital motion is calculated by solving a set of coupled ODEs numerically and then the waveforms can be calculated. With the wavefroms constructed we can then also determine the Signal to Noise ratio, using the noise models as described in [Robson 2018](https://arxiv.org/abs/1803.01944).

## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This code is written in FORTRAN with a [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. **Other compilers have not been tested.** The gfortran installation binaries can be found [here](https://gcc.gnu.org/wiki/GFortranBinariels), although typically gfortran comes pre-installed on most Linux/Unix systems. If you have [Homebew](https://brew.sh/) installed on OSX, you can simply run 


```
brew install gcc
```



### Starting steps
After [cloning the repo](https://help.github.com/en/articles/cloning-a-repository), the first thing to do is to set the path to the output files that the code will produce.This can be found in `src/parameters.f` 


```
echo 'export SCDir="/Users/tomkimpson/Data/SpinCurv/"' >> ~/.bash_profile
source ~/.bash_profile
```
Just change the path to some appropriate local destination

You can check the environemnt variable has been added to `bash_profile` by either `env` or `vim ~/.bashprofile`

Set this to point to a local direcory.

The code should then run as is, out of the box. Try

```
run.py
```

to compile and run the code. Once you have checked that everything is running OK, you can then start playing. The code structure (mdoules, subroutines etc.) is outlined below.


If making edits to the code, try to keep to the [FORTRAN Style Guide](https://www.fortran90.org/src/best-practices.html)

## Structure

### src/

`parameters.f` defines all the system parameters. That is, anything that needs changing (e.g. eccentricity, orbital period, BH mass) can be modified in this module


`constants.f` is for calculations with those parameters for use later in the code. It can effectively be ignored - no changes should be necessary to this file

`main.f` is where the code program is run from. After setting up the initial conditions (`initial_conditions.f`) it then goes on to integrate the equations and save the output (`rk.f`). 

In turn, `rk.f` calls `derivatives.f` and `tensors.f` to calculate e.g. the curvature tensors, ODEs and then integrates numerically.

### tools/

`plot_trajectory.py`. As on the tin. Can switch between 3d and 2d plotting.

`plot_ds.py`. Compares the spatial difference between the lamba=0 and lamba=1 cases. Useful for quickly eyeing Roemer delay. Requires interpolation.


A python wrapper has been provided to compile and run the code, `run.py`. We use a `-O3` optimization. [See the docs](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) for discussion on the optimization flags


## Numerical Method
We integrate the equations using a Runge-Kutta-Fehlberg algorithm with adaptive stepsize. See [Press et al.](https://dl.acm.org/citation.cfm?id=141273)





## Authors

* **Tom Kimpson** 





