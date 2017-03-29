Introduction
============
This is a highly efficient C implementation of [Yuval Tassa's iLQG algorithm](http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization) (see also his [paper](https://homes.cs.washington.edu/~todorov/papers/TassaICRA14.pdf)) plus some improvements and extentions.

Efficiency is achieved using the following features:
* Problem-specific code (system function, cost function and their derivatives) is generated from an analytic problem description in the [Maxima language](http://maxima.sourceforge.net/)
* Auxillari variables and their derivatives used in the system function and cost function are reused avoiding unneccesary recalculations
* For symmetric matrices only the upper triangle is calculated
* There is an experimental implementation of multi-threading execution (activated via the -DMULTI_THREADED=1 compiler switch). But that didn't improve speed on my machine. What is still missing is the implementation of a multi-threaded line-search.

On my computer, this yields a speed improvement of 8ms per iteration vs. 1500ms of the original MATLAB implementation for the car-parking example from Yuval's MATLAB code.

As improvements I have added the handling of state dependent input constraints and an experimental regularisation via a modified cholesky decomposition.

Everything is still very raw, partially tested and undocumented. But if you are interested I will be more than happy to supply addition information and help with problems. In this case, please just open new a github issue in my repo.

Code Generation
===============
The problem-specific code is generated from a problem definition in the [Maxima language](http://maxima.sourceforge.net/) using the maxima package `gentran`. Unfortunately, the gentran package as supplied by the current maxima release is not fully functional. Therefore you will have to download the package from my [my repo](https://github.com/jgeisler0303/maxima) and replace the files in the `maxima/share/contrib/gentran` folder of your maxima installation with the files from my repo. On my Linux computer I got gentran only to work using [Steel Bank Common Lisp (SBCL)](http://www.sbcl.org/) implementation. Unfortunately, this meant I had to recompile maxima from source. On Windows, SBCL is the default Lisp implementation supplied with the [wxMaxima](http://andrejv.github.io/wxmaxima/) installation.

Code generation and compilation is easy using the MATLAB/Octave function `make_iLQG.m`. This function takes as first argument the string name of the problem definition .mac-file and as optional second argument additional compiler switches. These can be used to get reasonable console output by setting the following defines: `-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1`. Further, the define `-DFULL_DDP=0` can be used to switch off the use of second order derivates of the system function during the backward-pass of the iLQG algorithm. Using the define `-DMULTI_THREADED=1` experimental multi-threaded calculation of derivates and backpass can be enabled.

Problem Description
-------------------
Your optimal control problem has to be defined in a maxima batch file, i.e. a text file with `.mac` suffix containing maxima expressions. In this file at least five variables plus an optional sixth have to be defined:
* `x`: a list of symbols for the state names (e.g. `x: [state1, state2, state3]`)
* `u`: a list of symbols for the system inputs (e.g. `x: [in1]`)
* `f`: the time discrete state transition function as an array index by the state name symbols. E.g.
```
    f[state1]: state2*ts;
    f[state2]: state3*ts+offset_parameter;
    f[state3]: in1*ts;
```
* `L`: the running cost function
* `F`: the final cost function (may not depend on inputs)
* `h`: an array of input contraint functions. For every `h[i]` it is assumed that `h[i]<0` and every constraint is enforced in order of ascending index of h. Each `h[i]` may only depend on one input with positive or negative one as coefficient.

The symbols for states and inputs may not be defined. In the expressions of `f`, `L`, `F` and `h` any undefined symbol appart from state and input symbols is considered to be a parameter. Any defined symbol is considered to be an auxillari value. Auxillari values and their derivatives are calculated separately and before their use in `f`, `L`, `F` or `h` or their derivatives. Auxillari values may depend on other auxillari values, states, inputs or parameters. If you use auxillaries, make sure to prefix them with an apostrophe. Otherwise, the auxillari definition will be substitued in to your expression and you end up without the benefits of auxillaries.

Getting the Example to Run
==========================
As an example the definition of the car parking problem from Yuval Tassa's iLQG implementation and the classic [brachistochrone problem](https://en.wikipedia.org/wiki/Brachistochrone_curve) can be found in the examples folder. These definitions can be compiled into an [Octave Mex function](https://www.gnu.org/software/octave/doc/interpreter/Getting-Started-with-Mex_002dFiles.html) or a MATLAB mex function using the `make_iLQG.m` MATLAB/Octave function with the following arguments: `make_iLQG('optDef{Brachi|Car}', '-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1 -DFULL_DDP=0')`. This is assuming you have Maxima, and for Octave the `mkoctfile` package installed and the `DDP-Generator` base path is in your MATLAB or Octave path. In MATLAB you may also have to `setenv('MAXIMA', 'path to maxima');`. Upon success, a subfolder named `{problem name}_gen_files` is created, containing the generated problem specific code and an empty `err.log` file. In the problem folder, the compiled `iLQG{prablem name}.mex` file is created. Now you can run the `testCar` or `testBrachi` demo script.
# DDP-generator-DriftCar
