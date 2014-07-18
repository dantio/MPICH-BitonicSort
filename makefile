NAME=tweetonic
MPICC?=mpicc -std=c99 -Wall -g
SRC=${NAME}.c
OUT=${NAME}.out
NP=4

all: ${NAME}

tweetonic: ${NAME}.c
	${MPICC} -o ${OUT} ${SRC}

clean:
	rm ${OUT}

check:
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./${OUT} "app"

run1:
	mpirun -np 1 ./${OUT} "app"
	
run2: 
	mpirun -np 2 ./${OUT} "app"
	
run4: 
	mpirun -np 4 ./${OUT} "app"
	
run8: 
	mpirun -np 8 ./${OUT} "app"
	
run16: 
	mpirun -np 16 ./${OUT} "app"
	
style:
	astyle --style=google --fill-empty-lines ./${SRC}

	
