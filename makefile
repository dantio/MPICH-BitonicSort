NAME=tweetonic
MPICC?=mpicc -std=c99 -Wall -g
SRC=${NAME}.c
OUT=${NAME}.out

all: ${NAME}

tweetonic: ${NAME}.c
	${MPICC} -o  ${OUT} ${SRC}

clean:
	rm ${OUT}

check:
	valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes ./${OUT} "app"

run:
	mpirun -np ${NP} ${HOSTFILE} ./${OUT} ${KEY} ${FILES}
	
style:
	astyle --style=google --fill-empty-lines ./${SRC}

	
