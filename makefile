NAME=tweetonic
MPICC?=mpicc -std=c99 -Wall
SRC=${NAME}.c Dictionary.c Tweet.c
OUT=${NAME}.out
NP=4

all: ${NAME}

tweetonic: ${NAME}.c
	${MPICC} -o ${OUT} ${SRC}

clean:
	rm ${OUT}

run: all
	mpirun -np ${NP} ./${OUT}
	
style:
	astyle --style=google --fill-empty-lines ./${SRC}
	
