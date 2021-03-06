################################################################################
#
# Makefile project only supported on Mac OS X and Linux Platforms)
#
################################################################################
include ./findcudalib.mk

# ---------------------------------------------------------------------
# Objects
# ---------------------------------------------------------------------
SOLVER_PATH=.
MAIN=cocos
H_FOLDER=host_core

HEADERS_CORE=atom.h\
    logic_variables.h\
	aminoacid.h\
	domain.h\
	protein.h\
	globals.h\
	input_data.h\
	mathematics.h\
	utilities.h\
	supervisor.h\
	mas_agent.h\
	structure_agent.h\
	coordinator_agent.h\
	worker_agent.h\
	search_engine.h\
	icm.h\
	constraint_store.h\
	constraint.h\
	propagator.h\
	distance.h\
	gibbs.h\
	montecarlo.h\
	cuda_energy.h\
	cuda_propagators.h\
	cuda_rmsd.h\
	#constraint_all_distant.h\
	energy.h\
	all_distant.h\
	sang.h\
	mang.h\
	k_angle_shuffle.h\
	k_rang.h\
	constraint_cg.h\
	constraint_fragment.h\
	constraint_shell.h\
	coordinator_agent.h\
	cuda_all_distant.h\
	cuda_cg.h\
	cuda_propagation.h\
	energy_seq.h\
	fragment.h\
	rmsd_fast.h\
	search.h\
	statistics.h\
	structure_agent.h\

HEADERS=$(HEADERS_CORE:%.h=${H_FOLDER}/%.h)

SOURCES=${MAIN}.cu $(HEADERS:%.h=%.cu) 
OBJECTS = $(SOURCES:%.cu=%.o)

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CC = g++
##CC = clang++
#NVCC = /usr/local/cuda/bin/nvcc 
NVCC = /usr/local/cuda-5.5/bin/nvcc

DEPEND = g++ -MM
CCOPT =    # -Wall -W -Wno-sign-compare -ansi -pedantic -Werror -Wno-deprecated
LINKOPT =

## debug, don't optimize 
CCOPT   += -m64 -DIL_STD

## debug, optimize
## use --ptxas-options=-v for cuda (verbose for each kernel)
## use -ftz=false -prec-div=true -prec-sqrt=true
## -L /usr/local/cuda/lib -lcurand
## –use_fast_math
CCOPT += -O3 
LINKOPT += -gx

PROFOPT = -pg -O2 -DNDEBUG

vpath %.o obj

## lib dirs -L...
CCLNDIRS= 
## include dirs -I...
INCLDIRS = -I$(SOLVER_PATH) -I$(SOLVER_PATH)/$(H_FOLDER) 

## Compiler Flags
CCFLAGS = $(CCOPT) $(INCLDIRS)
NVCCFLAGS = -arch=sm_35
#NVCCFLAGS = -arch=sm_21


all:	cocos

cocos: $(OBJECTS) 
	$(NVCC) $(NVCCFLAGS) $(CCLNDIRS) $(DFLAGS) -o $(MAIN) \
		$(OBJECTS:%=$(SOLVER_PATH)/obj/%)

$(OBJECTS): %.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $(CCFLAGS) $< -o $(SOLVER_PATH)/obj/$@

clean:
	rm -f obj/host_core/*.o 
	 rm -f obj/cocos.o
