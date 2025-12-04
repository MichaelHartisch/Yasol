# Yasol
Yasol is a search-based solver that is able to deal with multistage robust linear discrete optimization problems, with final mixed-integer recourse and a discrete uncertainty set, defined via linear constraints, which even can be decision-dependent.

# Building
Ideally, you only have to run

      ./configure [OPTIONS]
      cd build
      make check

See `./configure -h` for a full list of options.  
Depending on the selected configuration, the script may ask you to provide the path to your locally installed IP solver (CoinOR Tools, CPLEX, or HiGHS). If you do not have a solver installed, it can also download and install one for you (CoinOR Tools and HiGHS only).  
If you already have the solvers on your system, you may simply press **Enter** when prompted and let CMake attempt to locate the required installations automatically.  
The most common calls will look like this:  
Compiling with cplex: `./configure -c`

Compiling with highs: `./configure -g`

Compiling with coinor: `./configure -o`

After successfully configuring the solver, a `build` directory will appear in the project’s root folder.
Running `make` (or `make check` for a quick verification of the compiled program) will produce the yasol executable. The build process also generates the corresponding library files, which you can use if you wish to include Yasol in your own project.


# Running
From your `build` folder you can now run our solver by calling
      
      ./yasol [PATH_TO_INSTANCE] [OPTIONS]

You find some instances in the folder 'examples' in the main directory. A call then might look like this:

      ./yasol ../examples/knapsack.qlp --timeLimit=4 --isSimplyRestricted=1

# Yasol as Library
In case you want to use the solver as a callable library, you find for your convinence the folder myCMakeProject in the main directory. This folder contains a small project calling yasol. Depending on the IP solver you have compiled yasol with, you can now build this project as follows:

cplex: `cmake -B build -DUSE_CPLEX=ON -DYASOL_PATH=/path/to/yasol/main/directory`  
highs: `cmake -B build -DUSE_HIGHS=ON -DYASOL_PATH=/path/to/yasol/main/directory`  
cplex: `cmake -B build -DUSE_COIN=ON -DYASOL_PATH=/path/to/yasol/main/directory`  

Here, `YASOL_PATH` should be the main directory, **not** the `build` folder. Then do `cd build && make`.

The program `MyYasol' then builds the famous bilevel problem from Moore & Bard as QIP using our datastructures, calls the solver and outputs the result. 
# References

1. Hartisch, M. and Lorenz, U., 2022. A general model-and-run solver for multistage robust discrete linear optimization. arXiv preprint arXiv:2210.11132. [DOI](https://doi.org/10.48550/arXiv.2210.11132)
3. Hartisch, M., 2020. Quantified integer programming with polyhedral and decision-dependent uncertainty. PhD thesis. [DOI](https://doi.org/10.25819/ubsi/4841)
2. Ederer, T., Hartisch, M., Lorenz, U., Opfer, T. and Wolf, J., 2017. Yasol: an open source solver for quantified mixed integer programs. In Advances in Computer Games (pp. 224-233).  [DOI]( https://doi.org/10.1007/978-3-319-71649-7_19)


# Other Links
For further information please visit
[GitHub Pages](https://yasolqipsolver.github.io/yasol.github.io/)
or our [Homepage](http://www.q-mip.org/).
