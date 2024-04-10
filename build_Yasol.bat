#Make the changes according to your system here
#----------------------------------------------

#1. If neccessary change the compiler and archiver.
export YASOL_CC=/usr/bin/clang++
export YASOL_AR=/usr/bin/ar

#2. Depending on you system (Mac or Linux) ensure the correct lines are commented. They should be corrected after downloading the osx or linux version from www.q-mip.org
#If you are a Mac user the following lines should not be commented out
export YASOL_SYSTEM="x86-64_osx"
export YASOL_LIBFORMAT=""
export YASOL_OS="-DLINUX"

#If you are a Linux user the following lines should not be commented out
#export YASOL_SYSTEM="x86-64_linux"
#export YASOL_LIBFORMAT="-static-libgcc"
#export YASOL_OS="-DLINUX"

#If you are a Windows user: remove the comments here. 
#export YASOL_OS="-DWINDOWS"
#export YASOL_SYSTEM="x64_win64"


#3. Tell the compiler where to find the solver you want to use as external LP solver. Remove the comments at the appropriate lines. 
export YASOL_CPLEX_PATH=/Applications/CPLEX_Studio221/cplex
#export YASOL_CLP_PATH=/nethome/yasol/CbcClp

#There should be nothing you have to worry about below
#-----------------------------------------------------

OPTION_CPLEX=0
OPTION_CLP=0
OPTION_CGL=0
OPTION_DEBUG=0
if [[ -z "$@" ]]; then
    OPTION_CPLEX=1
    OPTION_CLP=0
    OPTION_CGL=0
    OPTION_DEBUG=0
    echo "Detected no additional options. Use standard setting, i.e. cplex, without CGL, optimized"
fi

for option in "$@"
do
    if [[ $option == "clp" ]]; then
        OPTION_CLP=1
        echo "Accepted Option: Use CLP as LP Solver"
    elif [[ $option == "cplex" ]]; then
        OPTION_CPLEX=1
        echo "Accepted Option: Use CPLEX as LP Solver"
    elif [[ $option == "cgl" ]]; then
        OPTION_CGL=1
        echo "Accepted Option: Use CGL as additional cut generator"
    elif [[ $option == "debug" ]]; then
        OPTION_DEBUG=1
        echo "Accepted Option: Compile in debug mode"
    else
        echo "Unknown Option: $option"
        exit
    fi
done

if [[ $OPTION_CLP == "1" ]] && [[ $OPTION_CGL == "1" ]]; then
    OPTION_CGL=0
    echo "Info: Removed cgl option, as using CLP as LP solver already includes the necessary library"
fi

if [[ $OPTION_CLP == "1" ]] && [[ $OPTION_CPLEX == "1" ]]; then
    echo "Error: Cannot use more than one LP solver"
    exit
fi

if [[ $OPTION_CLP == "0" ]] && [[ $OPTION_CPLEX == "0" ]]; then
    echo "Error: No LP solver specified. Try to use cplex as LP solver..."
    OPTION_CPLEX=1
fi

#Use Cplex
if [[ $OPTION_CPLEX == "1" ]]; then
    if [[ -z "$YASOL_CPLEX_PATH" ]]; then
        echo "Error: Set YASOL_CPLEX_PATH in build_Yasol.bat to the path of your CPLEX distribution."
        exit
    fi
    export YASOL_USED_SOLVER="-DCOMPILE_WITH_CPLEX_C -DUSE_CPLEX_C -DUSE_NBD_CPLEX_C"
    export YASOL_PATHS_SOLVER="-I"$YASOL_CPLEX_PATH"/include"
    export YASOL_LIB_SOLVER=$YASOL_CPLEX_PATH"/lib/"$YASOL_SYSTEM"/static_pic/libcplex.a"
    export YASOL_USE_CGL="-DNO_CGL"
fi


#Turn On CGL for external cuts from CGL at root (needs Coin OR)
if [[ $OPTION_CGL == "1" ]]; then
    if [[ -z "$YASOL_CLP_PATH" ]]; then
        echo "Error: Set YASOL_CLP_PATH in build_Yasol.bat to the path of your COIN-OR distribution when using the cgl option."
        exit
    fi
    export YASOL_USE_CGL=" "
    export YASOL_PATHS_CGL="-I"$YASOL_CLP_PATH"/include/coin -I"$YASOL_CLP_PATH"/dist/include/coin-or"
    export YASOL_LIB_CGL="-L"$YASOL_CLP_PATH"/dist/lib -L"$YASOL_CLP_PATH"/lib -lCbc -lOsiCbc -lCgl -lClp -lCoinUtils -lOsiClp -lOsi"
    export YASOL_OBJECT_CGL="src/ExternSolvers/QpExtSolCBC.o"
fi

#Use CLP Solver
if [[ $OPTION_CLP == "1" ]]; then
    if [[ -z "$YASOL_CLP_PATH" ]]; then
        echo "Error: Set YASOL_CLP_PATH in build_Yasol.bat to the path of your COIN OR distribution"
        exit
    fi
    export YASOL_USED_SOLVER="-DCOMPILE_WITH_CLP -DCOMPILE_WITH_CBC -DUSE_NBD_CLP -DUSE_CBC"
    export YASOL_PATHS_SOLVER="-I"$YASOL_CLP_PATH" -I"$YASOL_CLP_PATH"/dist/include/coin-or -I"$YASOL_CLP_PATH"/include/coin"
    export YASOL_LIB_SOLVER="-L"$YASOL_CLP_PATH"/dist/lib -L"$YASOL_CLP_PATH"/lib -lCbc -lOsiCbc -lCgl -lClp -lCoinUtils -lOsiClp -lOsi"
    export YASOL_OBJECT_SOLVER="src/ExternSolvers/QpExtSolCLP.o src/ExternSolvers/QpExtSolCBC.o"
fi

#Compiler/Debug Option
if [[ $OPTION_DEBUG == "1" ]]; then
    export YASOL_COMPILER_OPTION="-O1 -g -fsanitize=address -fno-omit-frame-pointer -fno-builtin-malloc"
else
    export YASOL_COMPILER_OPTION="-O3 -g"
fi

cd Solver
rm *.a
rm *.o
make -f makeSolver
cd ..
cd QBPSolver
rm *.a
rm *.o
make -f makeQBPSolver
cd ..
cd Yasol
rm *.o
make -f makeYasol
mv Yasol Debug/Yasol
ls -alt ../QBPSolver/*.a ../Solver/*.a Debug/Yasol*
#cd ..
