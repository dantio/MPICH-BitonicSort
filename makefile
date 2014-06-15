EXECS=tweetonic
MPICC?=mpicc -std=c99 -Wall
NP=4

all: ${EXECS}

tweetonic: ${EXECS}.c
	${MPICC} -o ${EXECS} ${EXECS}.c

clean:
	rm ${EXECS}

run:
	mpirun -np ${NP} ./${EXECS}
