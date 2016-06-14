# NOTE:
# This file is mainly useful when doing a manual installation of Armadillo.
# It is overwritten by CMake when doing an automatic installation.
#
# You may need to edit this file to reflect the type and capabilities
# of your system.  The defaults are for a Linux system and will need
# to be changed for other systems (e.g. Mac OS X).


CXX=g++
#CXX=mpicc
#CXX=g++-4.2
## Under MacOS you may have an old compiler as default (e.g. GCC 4.0).
## However, GCC 4.2 or later is available and preferable due to better
## handling of template code.

#CXX=CC
## When using the Sun Studio compiler

LIB_PATH= ./lib/
INCLUDE_PATH= ./include/
SRC_PATH= ./src/
BUILD_PATH= ./

##Environments
# Home unix
#ARMA_LIB_UNIX=-L/opt/OpenBLAS/lib -L/usr/lib/lapack -L/usr/share/doc -lm -llapack -lopenblas -larmadillo
#ARMA_LIB_UNIX=-L/usr/lib/atlas-base/atlas -L/usr/lib/lapack -L/usr/share/doc -llapack_atlas -lcblas -lm -llapack -lblas -larmadillo
#ARMA_LIB_UNIX=-L/usr/lib/libblas -L/usr/lib/lapack -L/usr/share/doc -lm -llapack -lblas -larmadillo
#ARMA_INCLUDE_FLAG_UNIX = -I/usr/include #-L/usr/lib/libblas -L/usr/lib/lapack -L/usr/share/doc
# New version of Armadillo. Keep the old versions in the repo to test against.
#ARMA_LIB_UNIX_STATIC=-L/usr/lib -L/usr/share/doc -L/usr/lib/atlas-base -L/usr/lib/atlas-base/atlas -larmadillo_static -lm -llapack -llapack_atlas -lcblas -latlas
#GLPK_LIB_FLAG_UNIX=-L/usr/local/lib/libglpk.so.35 -lglpk -lm
#GLPK_INCLUDE_FLAG_UNIX = -I$(LIB_PATH)glpk-4.52/src

# Windows MinGW
ARMA_INCLUDE_FLAG_WIN = -I$(LIB_PATH)armadillo-3.930.1/include/
ARMA_LIB_WIN=-L $(LIB_PATH)lib_win64/ -lblas -llapack -L$(LIB_PATH)armadillo-3.930.1/ -larmadillo.dll
# library dlls sourced from http://icl.cs.utk.edu/lapack-for-windows/lapack/#libraries_mingw
ARMA_LIB_WIN32=-L $(LIB_PATH)lib_win32/ -lblas -llapack -L$(LIB_PATH)armadillo-3.930.1_32/ -larmadillo.dll
#ARMA_INCLUDE_FLAG = -I $(LIB_PATH)armadillo-3.6.0/include
SHAPELIB_LIB_WIN=-L $(LIB_PATH)shapelib-1.3.0/ -lshp
SHAPELIB_INCLUDE_WIN=-I$(LIB_PATH)shapelib-1.3.0/
#GLPK_LIB_FLAG_WIN=-L$(LIB_PATH)winglpk-4.52/glpk-4.52/w64 -lglpk_4_52
#GLPK_INCLUDE_FLAG_WIN = -I$(LIB_PATH)winglpk-4.52/glpk-4.52/src
GLPK_LIB_FLAG_WIN=-L$(LIB_PATH)glpk_w64/lib -lglpk
GLPK_INCLUDE_FLAG_WIN = -I$(LIB_PATH)glpk_w64/include

# Cluster Red Hat Unix
#
ARMA_LIB_UNIX=-L/lib64/ -lm -llapack -lblas
ARMA_INCLUDE_FLAG_UNIX = -I$(LIB_PATH)armadillo/include/ -L/lib64/
ARMA_LIB_UNIX_STATIC=-L/usr/lib -L/usr/share/doc -L/usr/lib/atlas-base -L/usr/lib/atlas-base/atlas -larmadillo_static -lm -llapack -lblas -latlas
SHAPELIB_LIB_UNIX=-L $(LIB_PATH)shapelib-1.3.0/ -lshp
SHAPELIB_INCLUDE_UNIX=-I$(LIB_PATH)shapelib-1.3.0/
#GLPK_LIB_FLAG=-lglpk -lm
GLPK_LIB_FLAG_UNIX=-L$(LIB_PATH)glpk-4.55/lib/ -lglpk -lm
GLPK_INCLUDE_FLAG_UNIX = -I$(LIB_PATH)glpk-4.55/src
XSD_INCLUDE_UNIX= -I$(LIB_PATH)libxsd/ -I$(LIB_PATH)xerces-c-3.1.3/src/
XSD_LIB_UNIX= -L$(LIB_PATH)xerces-c-3.1.3/libxerces/lib/ -lxerces-c-3.1


PUGI_SRC_UNIX= $(LIB_PATH)pugixml/src/


UNIX_INCLUDE_FLAGS = -I$(INCLUDE_PATH) $(GLPK_INCLUDE_FLAG_UNIX) $(ARMA_INCLUDE_FLAG_UNIX) $(SHAPELIB_INCLUDE_UNIX) $(XSD_INCLUDE_UNIX) -I$(PUGI_SRC_UNIX)
WIN_INCLUDE_FLAGS = -I$(INCLUDE_PATH) $(ARMA_INCLUDE_FLAG_WIN) $(GLPK_INCLUDE_FLAG_WIN) $(SHAPELIB_INCLUDE_WIN) 

LIBSOURCES= $(SRC_PATH)smallmath.cpp $(SRC_PATH)ComputedObservable.cpp $(SRC_PATH)ParameterGroup.cpp  $(SRC_PATH)DnaMeasurement.cpp $(SRC_PATH)Residual.cpp $(SRC_PATH)GPSComputedObservable.cpp $(SRC_PATH)GPSBaseline.cpp  $(SRC_PATH)Station.cpp $(SRC_PATH)GPSResidual.cpp $(SRC_PATH)SearchQueue.cpp  $(SRC_PATH)MeasCycle.cpp  $(SRC_PATH)timer.cpp $(SRC_PATH)edge.cpp $(SRC_PATH)MeasSegment.cpp  $(SRC_PATH)MeasNetwork.cpp  $(SRC_PATH)L1Solver.cpp $(SRC_PATH)L1GLPKIPSolver.cpp  $(SRC_PATH)main.cpp $(SRC_PATH)L1GLPKSimplexDualSolver.cpp $(SRC_PATH)L1GLPKIPDualSolver.cpp $(SRC_PATH)L1GLPKSimplexSolver.cpp $(SRC_PATH)YClusterComputedObservable.cpp $(SRC_PATH)YClusterResidual.cpp $(SRC_PATH)YCluster.cpp $(SRC_PATH)FPointComputedObservable.cpp $(SRC_PATH)FPointResidual.cpp $(SRC_PATH)FPoint.cpp \
$(SRC_PATH)L2CGArmaSolver.cpp \
$(SRC_PATH)L2ArmaSolver.cpp \
$(SRC_PATH)DynaML-pimpl.cxx \
$(SRC_PATH)DynaML-pskel.cxx \
$(PUGI_SRC_UNIX)pugixml.cpp
# $(SRC_PATH)L1SimplexSolver.cpp $(SRC_PATH)L1ConvexSolver.cpp 

BOOST_INCLUDE_FLAG = -I /usr/include
## If you have Boost libraries, change /usr/include
## to point to where they are installed and uncomment
## the above line.


#EXTRA_LIB_FLAGS = -larmadillo #-llapck -lblas
## The above line is an example of the extra libraries
## that can be used by Armadillo. You will also need
## to modify "include/armadillo_bits/config.hpp"
## to indicate which extra libraries are present.
## If you're using Mac OS, comment out the line and
## instead use the line below.


#EXTRA_LIB_FLAGS = -framework Accelerate
## Uncomment the above line when using Mac OS
## and modify "include/armadillo_bits/config.hpp"
## to indicate that LAPACK and BLAS libraries
## are present


#EXTRA_LIB_FLAGS = -library=sunperf
## When using the Sun Studio compiler


LIB_FLAGS_UNIX = $(ARMA_LIB_UNIX) $(GLPK_LIB_FLAG_UNIX) $(SHAPELIB_LIB_UNIX) $(XSD_LIB_UNIX)
LIB_FLAGS_UNIX_STATIC = $(ARMA_LIB_UNIX_STATIC) $(GLPK_LIB_FLAG_UNIX) $(SHAPELIB_LIB_UNIX)
LIB_FLAGS_WIN = $(GLPK_LIB_FLAG_WIN) $(REGEX_LIB) $(ARMA_LIB_WIN) $(SHAPELIB_LIB_WIN) #$(GPSTK_LIB_FLAG_WIN) 
## NOTE: on Ubuntu and Debian based systems you may need to add 
## -lgfortran to LIB_FLAGS


OPT = -O2  -Wall -w -g
## TDM-specific : -static-libstdc++ is not a valid option.
## Apparently libstdc++ is statically linked by default
## 
## As the Armadillo library uses recursive templates,
## compilation times depend on the level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: good trade-off between compilation time and execution speed
## -O2: produces programs which have almost all possible speedups,
##      but compilation takes longer


#OPT = -xO4 -xannotate=no
## When using the Sun Studio compiler


#EXTRA_OPT = -fwhole-program
## Uncomment the above line if you're compiling 
## all source files into one program in a single hit.


#DEBUG = -DARMA_EXTRA_DEBUG
## Uncomment the above line to enable low-level
## debugging.  Lots of debugging information will
## be printed when a compiled program is run.
## Please enable this option when reporting bugs.


#FINAL = -DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## DANGEROUS!  Not recommended unless your code has been
## thoroughly tested.


#
#
#

UNIXCXXFLAGS = -g $(UNIX_INCLUDE_FLAGS) $(BOOST_INCLUDE_FLAG) $(DEBUG) $(FINAL) -std=c++11

WINCXXFLAGS = $(WIN_INCLUDE_FLAGS) $(DEBUG) $(FINAL) 

unix: $(BUILD_PATH)l1solve_64

static: $(BUILD_PATH)l1solve_64s

win: $(BUILD_PATH)l1solve_win64

win32: $(BUILD_PATH)l1solve_win32

			
$(BUILD_PATH)l1solve_64:
	$(CXX) $(UNIXCXXFLAGS) $(LIBSOURCES) \
	-m64 -o $@  $< -m64 $(OPT) $(EXTRA_OPT) $(LIB_FLAGS_UNIX)
			

$(BUILD_PATH)l1solve_64s:
	$(CXX) $(UNIXCXXFLAGS) $(LIBSOURCES) \
	-m64 -o $@  $< -m64 $(OPT) -static $(EXTRA_OPT) $(LIB_FLAGS_UNIX_STATIC)
			
			
$(BUILD_PATH)l1solve_win64:
	$(CXX) -m64 $(WINCXXFLAGS) $(LIBSOURCES) \
	-m64 -o $@  $< -m64 $(OPT) $(EXTRA_OPT) $(LIB_FLAGS_WIN)
			

$(BUILD_PATH)l1solve_win32:
	$(CXX) $(WINCXXFLAGS) $(LIBSOURCES) \
	-o $@  $< $(OPT) $(EXTRA_OPT) $(LIB_FLAGS_WIN)
			

.PHONY: clean

clean:
	rm -f $(BUILD_PATH)l1solve_64
	rm -f $(BUILD_PATH)l1solve_64s
	rm -f $(BUILD_PATH)l1solve_win64.exe
	rm -f $(BUILD_PATH)l1solve_win32.exe
