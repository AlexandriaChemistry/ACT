CC = clang++
LD = clang++

CFLAGS = -g -Wall -std=c++11
LDFLAGS =

RM = /bin/rm -f

OBJS = Crossover.o FitnessComputer.o GeneticAlgorithm.o ga_helpers.o Initializer.o Mutator.o ProbabilityComputer.o Selector.o Sorter.o Terminator.o test.o
GEN_HFILES = Initializer.h FitnessComputer.h ProbabilityComputer.h Crossover.h Mutator.h Terminator.h
EXEC = test

all: $(EXEC)

$(EXEC): $(OBJS)
	$(LD) -o $(EXEC) $(OBJS) $(LDFLAGS)

Crossover.o: Crossover.cpp Crossover.h aliases.h
	$(CC) $(CFLAGS) -c Crossover.cpp

FitnessComputer.o: FitnessComputer.cpp FitnessComputer.h aliases.h
	$(CC) $(CFLAGS) -c FitnessComputer.cpp

GeneticAlgorithm.o: GeneticAlgorithm.cpp GeneticAlgorithm.h aliases.h ga_helpers.h $(GEN_HFILES)
	$(CC) $(CFLAGS) -c GeneticAlgorithm.cpp

ga_helpers.o: ga_helpers.cpp ga_helpers.h aliases.h
	$(CC) $(CFLAGS) -c ga_helpers.cpp

Initializer.o: Initializer.cpp Initializer.h aliases.h
	$(CC) $(CFLAGS) -c Initializer.cpp

Mutator.o: Mutator.cpp Mutator.h
	$(CC) $(CFLAGS) -c Mutator.cpp

ProbabilityComputer.o: ProbabilityComputer.cpp ProbabilityComputer.h aliases.h
	$(CC) $(CFLAGS) -c ProbabilityComputer.cpp

Selector.o: Selector.cpp Selector.h aliases.h
	$(CC) $(CFLAGS) -c Selector.cpp

Sorter.o: Sorter.cpp Sorter.h aliases.h ga_helpers.h
	$(CC) $(CFLAGS) -c Sorter.cpp

Terminator.o: Terminator.cpp Terminator.h aliases.h ga_helpers.h
	$(CC) $(CFLAGS) -c Terminator.cpp

test.o: test.cpp aliases.h ga_includes.h
	$(CC) $(CFLAGS) -c test.cpp

clean:
	$(RM) $(EXEC) $(OBJS)