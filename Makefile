#assumes   module load comp/intel-15.0.3.187 
GSLROOT = /usr/local/other/SLES11.1/gsl/1.16/gnu-4.8.1
CXX = icpc
CC = icc

INCDIR = /opt/local/include $(GSLROOT)/include
LIBDIR = /opt/local/lib $(GSLROOT)/lib

LIBS  = gsl gslcblas m
CCFLAGS = -g -Wall -O3 -std=gnu99 

CC = gcc

OBJS = Subroutines.o TimePhaseMaximization.o LISAPathfinder.o

all: $(OBJS) mcmc test

LISAPathfinder.o: LISAPathfinder.c LISAPathfinder.h
	$(CC) $(CCFLAGS) -c LISAPathfinder.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

TimePhaseMaximization.o: TimePhaseMaximization.c TimePhaseMaximization.h
	$(CC) $(CCFLAGS) -c TimePhaseMaximization.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

Subroutines.o: Subroutines.c Subroutines.h 
	$(CC) $(CCFLAGS) -c Subroutines.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

mcmc: mcmc.c $(OBJS)
	$(CC) $(CCFLAGS) -o mcmc mcmc.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

test: test.c $(OBJS)
	$(CC) $(CCFLAGS) -o test test.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

clean:
	rm *.o mcmc test
