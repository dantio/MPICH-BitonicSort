MPICH-BitonicSort
=================

This is an implemantation of BitonicSort with MPICH.

Challange:
https://www.dropbox.com/s/5pbl1fs20rji80d/Belegarbeit.pdf

Parallel Solution by:
http://www.cs.rutgers.edu/~venugopa/parallel_summer2012/Resources/MPI_Bitonic_Graphic.png

HELP Operator
http://home.fhtw-berlin.de/~junghans/cref/CONCEPT/bit_shift.html

How to Run:
make
make run NP=N FILES=F KEY=K HOSTFILE="-f hostfile"
N = number of process
F = number of Files
K = search key

Ergebnise von anderen:
1:          2047757
2:          10269690

12345677:   3326929
12345678:   7907755
12345679:   22400396

23999999:   12720801
24000000:   13380831

Qsort: sort_local.c "as":
1:          2047757
2:          10269690

12345677:   6160299
12345678:   12725828
12345679:   7930934

23999999:   17163420
24000000:   9649720

Bitonic-Sort: tweetonic.c "as": 4NP
1:          2047757
2:          10269690

12345677:   16777989
12345678:   12048376
12345679:   14427356

23999999:   17163420
24000000:   9649720
