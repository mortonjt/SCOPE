CC=g++
DEBUG = -O3
OPTS = ${DEBUG} -Wall -std=c++0x -lpthread
all: TEST SCOPE++

SCOPE++: build_ghmm.o SCOPE.o Viterbi.o
	${CC} -o SCOPE++ build_ghmm.o SCOPE.o Viterbi.o ${OPTS} 

TEST: test.cpp build_ghmm.o Viterbi.o 
	${CC} -o TEST test.cpp build_ghmm.o Viterbi.o  ${OPTS}

build_ghmm.o: build_ghmm.cpp build_ghmm.h
	${CC} ${OPTS} -c build_ghmm.cpp

SCOPE.o: SCOPE.cpp SCOPE.h build_ghmm.h Viterbi.h 
	${CC} ${OPTS} -c SCOPE.cpp

Viterbi.o: Viterbi.cpp Viterbi.h build_ghmm.h matrix.hh
	${CC} ${OPTS} -c Viterbi.cpp

###############
# The following three commands are for creating optimized executables
# with a name indicating OS platform.  Use only for creating executables
# for distribution.
# SCOPE++.centOS: build_ghmm.cpp SCOPA.cpp Viterbi.cpp
# 	${CC} -o SCOPE++.centOS build_ghmm.cpp SCOPA.cpp Viterbi.cpp -O3

# SCOPE++.osx: build_ghmm.cpp SCOPA.cpp Viterbi.cpp
# 	${CC} -o SCOPE++.osx build_ghmm.cpp SCOPA.cpp Viterbi.cpp -O3

# SCOPE++.ubuntu: build_ghmm.cpp SCOPA.cpp Viterbi.cpp
# 	${CC} -o SCOPE++.ubuntu build_ghmm.cpp SCOPA.cpp Viterbi.cpp -O3


###############
clean:
	rm *.o TEST SCOPE++
