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
	mpirun -np 2 ./${OUT}
	
run2: all
	mpirun -np 4 ./${OUT}
	
run4: all
	mpirun -np 8 ./${OUT}
	
run16: all
	mpirun -np 16 ./${OUT}
	
style:
	astyle --style=google --fill-empty-lines ./${SRC}
	
