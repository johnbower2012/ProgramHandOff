## output file name
MAIN=main
PARAMETERS=parameters

## suffixes used:
#.SUFFIXES: .cpp .o .x .h

## software directories
#### include: *.h  obj: *.o  lib: *.a  src: *.cpp
SOFT=../software
IDIR=$(SOFT)/include
ODIR=$(SOFT)/build
LDIR=$(SOFT)/lib
SDIR=$(SOFT)/src

## compiler and compiler flags
CC=g++
CFLAGS=-O3 -Wall -g -std=c++11
AFLAGS=-I$(IDIR) #-fast -W -Wall -WShadow -Wconversion

## additional libraries to link
LIBS=

## dependecies
#### pattern substition adds folder location to *.h files
#_DEPS=coshfunc.h analysis.h system.h emulator.h mcmc.h parametermap.h
ALL_DEPS=system.h parametermap.h balancemodel.h analysis.h emulator.h mcmc.h coshfunc.h NN_parts.h
ALLDEPS=$(patsubst %, $(IDIR)/%, $(MAIN_DEPS))

## object files
#### pattern substition adds folder location to *.o files
#_OBJ=coshfunc.o analysis.o system.o emulator.o mcmc.o parametermap.o $(MAIN).o
ALL_OBJ=system.o parametermap.o balancemodel.o analysis.o emulator.o NN_parts.o mcmc.o coshfunc.o $(MAIN).o
ALLOBJ=$(patsubst %, $(ODIR)/%, $(MAIN_OBJ))

## dependecies
#### pattern substition adds folder location to *.h files
#_DEPS=coshfunc.h analysis.h system.h emulator.h mcmc.h parametermap.h
MAIN_DEPS=system.h parametermap.h balancemodel.h analysis.h emulator.h mcmc.h coshfunc.h NN_parts.h
MAINDEPS=$(patsubst %, $(IDIR)/%, $(MAIN_DEPS))

## object files
#### pattern substition adds folder location to *.o files
#_OBJ=coshfunc.o analysis.o system.o emulator.o mcmc.o parametermap.o $(MAIN).o
MAIN_OBJ=system.o parametermap.o balancemodel.o analysis.o emulator.o NN_parts.o mcmc.o coshfunc.o $(MAIN).o
MAINOBJ=$(patsubst %, $(ODIR)/%, $(MAIN_OBJ))

## dependecies
#### pattern substition adds folder location to *.h files
PARAMETERS_DEPS=system.h parametermap.h coshfunc.h
PARAMETERSDEPS=$(patsubst %, $(IDIR)/%, $(MAIN_DEPS))

## object files
#### pattern substition adds folder location to *.o files
PARAMETERS_OBJ=system.o parametermap.o coshfunc.o $(PARAMETERS).o
PARAMETERSOBJ=$(patsubst %, $(ODIR)/%, $(MAIN_OBJ))

