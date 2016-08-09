
INCDIR = /opt/local/include
LIBDIR = /opt/local/lib

LIBS  = gsl gslcblas m
CCFLAGS = -O3 -Wall -std=gnu99

CC = gcc

OBJS = BayesLine.o Subroutines.o TimePhaseMaximization.o LISAPathfinder.o Spacecraft.o

all: $(OBJS) mcmc test

BayesLine.o : BayesLine.c BayesLine.h
	$(CC) $(CCFLAGS) -c BayesLine.c 

LISAPathfinder.o: LISAPathfinder.c LISAPathfinder.h
	$(CC) $(CCFLAGS) -c LISAPathfinder.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

Spacecraft.o: Spacecraft.c Spacecraft.h
	$(CC) $(CCFLAGS) -c Spacecraft.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

TimePhaseMaximization.o: TimePhaseMaximization.c TimePhaseMaximization.h
	$(CC) $(CCFLAGS) -c TimePhaseMaximization.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

Subroutines.o: Subroutines.c Subroutines.h Spacecraft.h 
	$(CC) $(CCFLAGS) -c Subroutines.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

mcmc: mcmc.c $(OBJS)
	$(CC) $(CCFLAGS) -o mcmc mcmc.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

test: test.c $(OBJS)
	$(CC) $(CCFLAGS) -o test test.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

clean:
	rm *.o mcmc test
