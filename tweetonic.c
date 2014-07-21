/*
 ============================================================================
 Name        : Tweetonic.c
 Author      : Daniil Tomilow s0531603
               Eric Wündisch
 ============================================================================
 */

#include <mpi.h>        // MPI Library
#include <stdio.h>      // Printf
#include <stdlib.h>     // calloc, malloc
#include <string.h>     // memcmp
#include <time.h>       // Timer

// Bitonic sort
#define ASCENDING  1
#define DESCENDING -1

// Tweets
#define MAX_LINE_SIZE 1000

#define TSIZE 32
#define TNUM 100 //24000000 // Zeilen

#define FIN  "/mpidata/parsys14/gross/twitter.data."
#define FOUT "twitter.sort."

char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

int FNUM;
int linesToRead;

double brutto_start, brutto_end, netto_start, netto_end;

// Global Tweets array
char **TWEETS;

/**
 * Print error and exit programm
 */
void handle_error(const char* msg) {
    MPI_Finalize();
    fprintf(stderr, "Error: %s\n", msg);
    exit(255);
}

/**
 * Write ordered tweet to file
 */
void writeOrderedTweets(int rank, char **tweet, int size) {
    char f[20];
    sprintf(f, FOUT"%d", rank);
    
    remove(f);
    FILE* fu = fopen(f, "a");
    for (int i = 0; i < size; i++) {
        fprintf(fu, "%d %d %d\n", (short) tweet[i][0], (int) tweet[i][2],  tweet[i][6]);
    }
    
    fclose(fu);
}

/**
 * Swap Tweet
 */
static inline void swap(int i, int j) {
    char c[TSIZE];
    memcpy(c, TWEETS[i], TSIZE);
    memcpy(TWEETS[i], TWEETS[j], TSIZE);
    memcpy(TWEETS[j], c, TSIZE);
}

/**
 * Compare Tweet
 */
int compare(const void* ptr1, const void* ptr2) {
	int i;
	char* t1 = *((char**) ptr1);
	char* t2 = *((char**) ptr2);
	for (i=6; i<TSIZE; i++) {
		if (t1[i] > t2[i]) return -1;
		if (t2[i] > t1[i]) return 1;
	}
	return 0;
}   

/**
 * Bitonic sort of non power of two
 */
void bitonicMerge(int lo, int n, int dir) {
    if ( n > 1) {
        int m = n / 2;
        for (int i=lo; i<lo+m; i++) if(dir == compare(&TWEETS[i], &TWEETS[i + m])) swap(i, i + m);
        bitonicMerge(lo, m, dir);
        bitonicMerge(lo + m, m, dir);
    }
}

void bitonicSort(int lo, int n, int dir) {
    if (n > 1) {
        int m = n / 2;
        bitonicSort(lo, m, ASCENDING);
        bitonicSort(lo + m, m, DESCENDING);
        bitonicMerge(lo, n, dir);
    }
}

/**
BITONIC WITH POWER OF TWO
void bitonic(unsigned int lines) {
    unsigned int i,j,ij,k,c;
    for (k = 2; k <= lines; k = 2 * k)
        for (j = k >> 1; j > 0; j = j >> 1)
            for (i = 0; i < lines; i++)
                if((ij = i ^ j) < lines)
                    if (ij > i && (c = compare(TWEETS[i], TWEETS[ij])))
                        if ((!(i & k) && c == 1) || ( (i & k) && c == -1))
                            swap(i,ij);
}
*/

/**
 * FROM FORTENBACHER
 */
unsigned int readNumber(char** lptr) {
    char* ptr = *lptr;
    char* line = *lptr;
    while (*ptr != ' ') ptr++;
    *ptr = 0;
    *lptr = ptr+1;
    return (unsigned int)atoi(line);
}

unsigned int readMonth(char** lptr) {
    char* ptr = *lptr;
    char* line = *lptr;
    while (*ptr != ' ') ptr++;
    *ptr = 0;
    *lptr = ptr+1;
    unsigned int i;
    unsigned int m;
    for (i=0, m=1; i<12; i++, m++)
        if (strncmp(line, MONTHS[i], 3) == 0)
            return m;
            
    handle_error("invalid month");
    return -1;
}

unsigned int countHits(const char* line, const char* key) {
    int n = strlen(key);
    int k = strlen(line) - n;
    int i;
    unsigned int hits = 0;
    for (i = 0; i < k; i++, line++) {
        if (*line == *key) {
            if (memcmp(line, key, n) == 0)
                hits++;
        }
    }
    return hits;
}

void writeTweet(char* tweet, const int fn, const int ln, const int hits,
		const int month, const int day, char* line) {
  short* ptr1 = (short*) tweet;
  *ptr1 = (short) fn;
  int* ptr2 = (int*) (tweet+2);
  *ptr2 = ln;
  *(tweet+6) = (char) hits;
  *(tweet+7) = (char) month;
  *(tweet+8) = (char) day;
  int i;
  int n = TSIZE-9;
  for (i=strlen(line); i<n; i++) line[i] = ' '; // padding
  memcpy(tweet+9, line, n);
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

/**
 * ALLOCALTE ONE BLOG MEMORY FOR SPECIFIED LINES
 */
char **allocTweets(int lines) {
    char **tweets = calloc(lines, sizeof(char*)); // Array of pointers
    char *data    = calloc(lines, sizeof(char) * TSIZE); // The data
    
    for (int i = 0; i < lines ; i++) tweets[i] = data + i * TSIZE;
    
    return tweets;
}

/**
 * COPY TWEET DATA FROM RIGHT TO LEFT
 */
static inline void copyTweet(char *t1, char *t2) {
    memcpy(t1, t2, TSIZE);   
}

/**
 * MAIN PARALLEL EXECUTION
 */
void parallel(FILE* files[], const char* key, int rank, int size, int readedLines) {
    // First we go down
    int down = 1;
    for(int next = size > 1; next > 0 && next < size ; ) {
    
        // Calculate who is sender and receiver
        // 4 Proc - A B C D E F G H
        // Next 1:
        // A -> B  C -> D  E -> F  G -> H
        // Next 2:
        // A -> C  B -> D  C - > E  D -> F  E -> G  F -> H
        // etc.
        int sender   = (next == 1 && (!isLastProc(rank, size) && (rank == 0 || rank % 2 == 0)))
                       || (next > 1 && (rank < size - next) );
                       
        int receiver = (next == 1 && ( isLastProc(rank, size) || (rank != 0 && rank % 2 != 0)))
                       || (next > 1 && rank > next -1 );
        ///////////////////////////////////////////////////
        // Compare HEIGH
        ///////////////////////////////////////////////////
        if(sender) {
            printf("1.RANK %d SEND LAST TWEET FROM %d\n",rank, rank + next);
            // Send last (worst) tweet to right
            MPI_Request r;
	        MPI_Send(TWEETS[linesToRead - 1], TSIZE, MPI_CHAR, rank + next, 0, MPI_COMM_WORLD);
            
            MPI_Status status;
            int  nbytes = 0;
            MPI_Probe(rank + next, 2, MPI_COMM_WORLD, &status);
            if(MPI_Get_count(&status, MPI_CHAR, &nbytes) == MPI_UNDEFINED) {
                handle_error("Error msg");
            }
            
            // kopie von nbytes, da es nach MPI_Recv komische zahl hat
            int tweets = nbytes / TSIZE;
            
            printf("5.RANK %d GET %d TWEET FROM %d\n",rank, tweets, rank + next);
  
            if(tweets != 0) {
                char **tweetsFromRight = allocTweets(tweets);
                
                MPI_Recv(tweetsFromRight, tweets * TSIZE, MPI_CHAR, rank + next, 2, MPI_COMM_WORLD, &status);
            
                int j = linesToRead - 1; // Iterator goes from bad to good tweets
                int g = 0;
                
                // Welche Tweets sind  besser?
                int i = 0;
                for(; j > -1 && g < tweets; j--, i++) {
                    if(compare(&tweetsFromRight[g], &TWEETS[j]) == -1) {
                        g++;
                    } else {
                        break;
                    }
                }
                

                
                printf("6.RANK %d SEND %d TWEET TO %d\n",rank, i, rank + next);
                MPI_Send(TWEETS, i * TSIZE, MPI_CHAR, rank + next, 3, MPI_COMM_WORLD);
                
                brutto_start = MPI_Wtime();
                // Build bitonic sequenz
                int h = linesToRead - i; // Iterator goes from bad tweets - tweets that are better
                for(int k = i - 1; k  > -1 ; k--, h++) {
                    copyTweet(TWEETS[h], tweetsFromRight[k]);
                }
                brutto_end += MPI_Wtime() - brutto_start;
                
                free(tweetsFromRight[0]);
                free(tweetsFromRight);
                
                brutto_start = MPI_Wtime();
                bitonicSort(0, readedLines, ASCENDING);
                brutto_end += MPI_Wtime() - brutto_start;
            }
        }
        ///////////////////////////////////////////////////
        // Compare Low
        ///////////////////////////////////////////////////
        if(receiver) {
            char getData[TSIZE];
            printf("2.RANK %d GET LAST TWEET FROM %d\n", rank, rank - next);
            MPI_Recv(getData, TSIZE, MPI_CHAR, rank - next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            int betterTweets = 0;
            for(int i = 0; betterTweets < linesToRead; i++) {
                int c = compare( TWEETS[i], getData );
                if( c == -1 )
                    betterTweets++;
                else
                    break;
            }
            printf("OK %d\n", betterTweets);
            if(betterTweets > 0) {
                MPI_Request r;
                MPI_Isend(TWEETS, betterTweets * TSIZE, MPI_CHAR, rank - next, 2, MPI_COMM_WORLD, &r);
                printf("3.RANK %d SEND %d TWEET TO %d\n",rank, betterTweets, rank - next);
                
                MPI_Status status;
                int  nbytes = 0;
                MPI_Probe(rank - next, 3, MPI_COMM_WORLD, &status);
                if(MPI_Get_count(&status, MPI_CHAR, &nbytes) == MPI_UNDEFINED) {
                    handle_error("ERROR: MPI_Get_count TAG 3");
                }
                
                int tweets = nbytes  / TSIZE;
                char **tweetsFromLeft = allocTweets(tweets);
                
                // Wir bekommenn die anderen Tweets zurueck
                MPI_Recv(getData, tweets * TSIZE, MPI_CHAR, rank - next, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("4.RANK %d GET %d TWEET FROM %d\n",rank, tweets, rank  - next);

                for(int k = 0; k < tweets; k++) {
                    copyTweet(TWEETS[k], tweetsFromLeft[k]);
                }
                
                free(tweetsFromLeft[0]);
                free(tweetsFromLeft);
                
            } else {
                // Ignore
                MPI_Request request;
                MPI_Isend(NULL, 0, MPI_INT, rank - next, 2, MPI_COMM_WORLD, &request);
            }
        }
        
        if(next > (size / (size / 2)) - 1 ) {
            // Go UP!
            down = 0;
        }
        
        if(down) next++;
        else next--;
        
	/*
        brutto_start = MPI_Wtime();
        bitonicSort(0, readedLines, ASCENDING);
        brutto_end += MPI_Wtime() - brutto_start;
	*/
    }  
}

/**
 * MAIN EXECUTURE
 */
void exec(const int numFiles, const int rank, const int size, const char* key) {
    
    FILE *files[FNUM];
    
    char fileName[20];
    linesToRead = TNUM / size;
    
    // Allocate size for tweets
    TWEETS = allocTweets(FNUM * linesToRead * size);
    char *willy = TWEETS[0];
    netto_start = MPI_Wtime();
    
    int iLine = 0;
    for(int f = 0; f < FNUM; f++) {
        sprintf(fileName, FIN"%d",f);
        files[f] = fopen(fileName, "rb");
        if (files[f] == NULL) handle_error("Cannon open File .\n");
        
        int globalstart = rank * linesToRead;
        int globalend   = globalstart + linesToRead;
        int sortStart = iLine;
        
        // calc lines to read for last proc
        if (isLastProc(rank, size)) {
            globalend = TNUM;
            linesToRead = TNUM - linesToRead * (size -1);
        };
        
        char buf[MAX_LINE_SIZE];
        char* line;
        
        unsigned int fn;
        unsigned int ln;
        unsigned int month;
        unsigned int day;
        unsigned int hits;

        for(int lines = 0; iLine < linesToRead * FNUM && (line = fgets(buf, MAX_LINE_SIZE, files[f])) != NULL; ++lines) {
            if(lines < globalstart) continue;
            if(lines > globalend -1) break;
            
            fn = readNumber(&line);
            
            ln = readNumber(&line);
            
            month = readMonth(&line);
            day = readNumber(&line);
            
            hits = countHits(line, key);
            
            writeTweet(TWEETS[iLine], fn, ln, hits, month, day, line);

            iLine++;
        }
	
        // sort tweets
        brutto_start = MPI_Wtime();
        qsort(TWEETS, iLine, sizeof(char*), compare);
        //bitonicSort(sortStart, iLine, ASCENDING);
        brutto_end = MPI_Wtime() - brutto_start;

	
        // wait till all proc are finished
        //MPI_Barrier(MPI_COMM_WORLD);
        
        // change tweets between procs
        //parallel(files, key, rank, size, iLine);
	if(rank == 0) printf("\nFile %d ready", f);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    //bitonicSort(0, linesToRead, ASCENDING);
    
    netto_end = MPI_Wtime() - netto_start;
    
    // Write tweets to file
    writeOrderedTweets(rank, TWEETS, linesToRead);
    
    free(willy);
    free(TWEETS);
    
    // Close files
    for(int i = 0; i < FNUM; i++) fclose(files[i]);
    
    if(rank != 0) {
        // Send Time to master
        double sendTime[2] = { netto_end, brutto_end };
        double buff;
        MPI_Gather( &sendTime, 2, MPI_DOUBLE, &buff, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    }
}

/**
 * PRINT TWEET BY SPECIFIED RANK
 */
void printTweetAtRank(int tweetRank, int size) {
    int platz;
    if(tweetRank < 1) {
        platz = TNUM + tweetRank + 1 ;
        tweetRank = TNUM / size - tweetRank + 1;
    } else {
        platz = tweetRank;
    }
    
    int f = 0;
    
    while(tweetRank > (TNUM / size)) {
        tweetRank -= TNUM / size;
        f++;
    }
    char fileName[20];
    sprintf(fileName, FOUT"%d",f);
    
    FILE *sortFile = fopen(fileName, "r");
    if(sortFile == NULL) handle_error("Error in printTweetAtRank\n");
    
    char buf[MAX_LINE_SIZE];
    char* line;
    
    // Read from sorted File
    for(int lines = 0; (line = fgets(buf, MAX_LINE_SIZE, sortFile)) != NULL; ++lines) {
        if(lines == tweetRank) {
            int fn = fn = readNumber(&line);
            int ln = readNumber(&line);
            
            sprintf(fileName, FIN"%d",fn);
            // Read from origin twitter file
            FILE *file = fopen(fileName, "r");
            if(file == NULL) handle_error("Error in printTweetAtRank\n");
            
            for(lines = 0; (line = fgets(buf, MAX_LINE_SIZE, file)) != NULL; lines++) {
                if(lines  == ln) {
                    printf("Tweet auf Platz: %d\n===\n", platz);
                    printf("%s ===\n", line);
                    //print("===\n");
                    fclose(file);
                    break;
                }
            }
            fclose(sortFile);
            break;
        }
    }
}

int main(int argc, char** argv) {

    if(argc != 3) {
        fprintf(stderr, "make run NP=X KEY=Y FILES=Z\n");
        return EXIT_FAILURE;
    }
    FNUM = atoi(argv[2]);
    
    // Initialize MPI
    MPI_Init(NULL, NULL);
    
    // Get the number of processes and rank
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0) {
        printf("Number of process: '%d'\nFiles: '%d'\nKey: '%s'\n=== Start execution", size, FNUM, argv[1]);
    }
    
    // Execute main programm
    exec(FNUM, rank, size, argv[1]);
    
    if(rank == 0) {
        printf("\n Zeitmessung:  \n");
        printf("Rank: %d  - Netto ='%f'; Brutto ='%f'\n",rank, netto_end, brutto_end);
        
        double getTime[size][2];
        double send;
        
        MPI_Gather( &send, 1, MPI_DOUBLE, getTime, 2, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
        
        for(int i = 1; i < size; i++) {
            printf("Rank: %d  - Netto ='%f'; Brutto ='%f'\n",i, getTime[i][0], getTime[i][1]);
        }
        printf("\n == Erfolgreich abgeshlossen == \n");
        printf("Das Suchwort war %s \n", argv[1]);
        int p = -1;
        
        while(1) {
            printf("Zeige Tweet mit Relevanz ( '0' to exit, '-1' Get Last ) : \n");
            scanf("%d", &p);
            if(p == 0) break;
            printTweetAtRank(p, size);
        }
    }
    
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}
