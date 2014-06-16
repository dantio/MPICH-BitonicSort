EXECS=tweetonic
MPICC?=mpicc -std=c99 -Wall
NP=4

all: ${EXECS}

tweetonic: ${EXECS}.c
	${MPICC} -o ${EXECS}.out ${EXECS}.c Dictionary.c

clean:
	rm ${EXECS}.out

run: all
	mpirun -np ${NP} ./${EXECS}.out
	
style:
	astyle --style=google --fill-empty-lines ./${EXECS}.c
	
