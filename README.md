# Yasol
Yasol is a search-based solver that is able to deal with multistage robust linear discrete optimization problems, with final mixed-integer recourse actions and a discrete uncertainty set, which even can be decision-dependent.

# Install and Run (on Mac and Linux)
Open the build_Yasol.bat to edit it.
1. If neccessary change the compiler and archiver.
2. Depending on you system (Mac or Linux) ensure the correct lines are (un)commented. 
3. Tell the compiler where to find the solver you want to use as external LP solver (clp/cbc and/or cplex). This might look like this:

export YASOL_CLP_PATH=/opt/tools/coinor_tools\
export YASOL_CPLEX_PATH=/opt/ibm/ILOG/CPLEX_Studio1261/cplex\
or\
export YASOL_CLP_PATH=/nethome/user/CbcClp\
export YASOL_CPLEX_PATH=/Applications/CPLEX_Studio221/cplex \

Start the script by typing "./build_Yasol.bat"
You can use the following options
cplex  - to use cplex as lp solver
clp    - to use clp as lp solver
cgl    - to add cgl for obtaining additional cuts; this options als need the YASOL_CLP_PATH
debug  - for having more information when having to debug. Not recommended for inexperienced users.

Calls might look like this:
./build_Yasol.bat cplex
./build_Yasol.bat clp
./build_Yasol.bat cplex cgl

If you are using CPLEX and are working on a mac with next gen chip we recommend preprending "arch -x86_64"

If the process ended without any error you can try to run Yasol. You will find the executable in the Yasol/Debug folder.
You can solve an instance by calling Yasol via ".\Yasol MyInstance.qlp"

When trying to execute Yasol you might encounter library errors. You may need to
  - add 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/PATH/TO/CbcClp/dist/lib' to your ~/.bashrc (Linux)
  - add 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/PATH/TO/CbcClp/dist/lib' to ~/.bashrc (OS X)

We also observed that we had to add /usr/bin/x86_64-linux-gnu to the LD_LIBRARY_PATH.

# Install and Run (Windows)
The latest version of Yasol is not supported on Windows. However, older version of our solver exist, that can run on Windows. Please see the instructions on the [GitHub Pages](https://yasolqipsolver.github.io/yasol.github.io/)
or our [Homepage](http://www.q-mip.org/). 
# Yasol.ini
The Yasol.ini file can be used to change some high level settings of the solver. We recommend that you change them only if you really know what is going to happen. Always feel free to reach out to us.
One special note regarding the Yasol.ini: The standard setting we selected is "isSimplyRestricted=0". This makes the search rather slow and careful when there is a universal constraint system. 
Setting isSimplyRestricted=1 might significantly speed up the search, but this is only allowed if the instances has special properties. If your instance does not have the demanded property this setting will
result in wrong results. Please see our publications or feel free to contact us any time.

# Other Links
For further information please visit
[GitHub Pages](https://yasolqipsolver.github.io/yasol.github.io/)
or our [Homepage](http://www.q-mip.org/).
