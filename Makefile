
INCDIR = #/opt/local/include
LIBDIR = #/opt/local/lib

LIBS  = gsl gslcblas m
CCFLAGS = -g -Wall -O3 -std=gnu99

CC = gcc

OBJS = Subroutines.o TimePhaseMaximization.o

all: $(OBJS) mcmc


TimePhaseMaximization.o: TimePhaseMaximization.c TimePhaseMaximization.h
	$(CC) $(CCFLAGS) -c TimePhaseMaximization.c

Subroutines.o: Subroutines.c Subroutines.h
	$(CC) $(CCFLAGS) -c Subroutines.c 

mcmc: mcmc.c 
	$(CC) $(CCFLAGS) -o mcmc mcmc.c $(OBJS) $(LIBS:%=-l%)

clean:
	rm *.o mcmc
