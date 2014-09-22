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

//#define DEBUG 1

// Bitonic sort
#define ASCENDING  1
#define DESCENDING -1

// TAGS
#define WORSTTWEET 0
#define NUMBEST 1
#define THREE 3
#define FOUR 4
#define FIVE 5

// Tweets
#define MAX_LINE_SIZE 1000

#define TSIZE 32
#define TNUM 240000 // 24000000 // Zeilen PRO FILE

//#define FIN  "/mpidata/parsys14/gross/twitter.data."
#define FIN  "/mpidata/parsys14/gross/twitter.data."
#define FOUT "twitter.sort."

#define S_INT 8

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

void printTweet(const char* t) {
    printf("[");
    int i;
    int* ptr = (int*) (t + 2);
    //printf("%d", *ptr);
    for (i=0; i<TSIZE; i++) {
        int k = t[i];
        printf("%2x ", k<0?k+256:k);
    }
    printf("]\n");
}

/**
 * Write ordered tweet to file
 */
void writeOrderedTweets(int rank, char **T, int size) {
    char f[20];
    sprintf(f, FOUT"%d", rank);
    char *tweet;
    remove(f);
    FILE* fu = fopen(f, "a");
    if(rank == 0 || rank %2 == 0) {
        for (int i = 0; i < size; i++) {
            tweet = T[i];
            fprintf(fu, "%d %d %d\n", *((short*)tweet), *((int*)(tweet + 2)), *(tweet +6));
        }
    } else {
        for (int i = size - 1; i > -1; i--) {
            tweet = T[i];
            fprintf(fu, "%d %d %d\n", *((short*)tweet), *((int*)(tweet + 2)), *(tweet +6));
        }
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
 * Compare Tweet ACSENDING
 */
static inline int compare(const void* ptr1, const void* ptr2) {
    char* t1 = (char*) ptr1;
    char* t2 = (char*) ptr2;
    
    for (int i = 6; i < TSIZE - 6; i++) {
        if (*(t1 + i) > *(t2 + i)) return -1;
        if (*(t2 + i) > *(t1 + i)) return 1;
    }
    
    return 0;
}

/**
 * Compare Tweet DESCENDING
 */
static inline int compareD(const void* ptr1, const void* ptr2) {
    char* t1 = (char*) ptr1;
    char* t2 = (char*) ptr2;
    
    for (int i = 6; i < TSIZE - 6; i++) {
        if (*(t1 + i) > *(t2 + i)) return   1;
        if (*(t2 + i) > *(t1 + i)) return  -1;
    }
    
    return 0;
}

static inline int greatestP(int n) {
    int k = 1;
    while(k < n) k <<= 1;
    return k >> 1;
}

/**
 * Bitonic sort of non power of two
 */
void bitonicMerge(int lo, int n, int dir) {
    if (n > 1) {
        int m = greatestP(n);
        for (int i = lo; i < lo + n - m; i++)
            if (dir == compare(TWEETS[i], TWEETS[i + m]))
                swap(i, i + m);
                
        bitonicMerge(lo, m, dir);
        bitonicMerge(lo + m, n - m, dir);
    }
}

void bitonicSort(int lo, int n, int dir) {
    if (n > 1) {
        int m = n / 2;
        bitonicSort(lo, m, (dir == 1 ? -1 : 1));
        bitonicSort(lo + m, n - m, dir);
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
                    if (ij > i && (c = compare(&TWEETS[i], &TWEETS[ij])))
                        if ((!(i & k) && c == 1) || ( (i & k) && c == -1))
                            swap(i,ij);
}
*/

/**
 * FROM FORTENBACHER
 */
int
readNumber(char** lptr) {
    char* ptr = *lptr;
    char* line = *lptr;
    while (*ptr != ' ') ptr++;
    *ptr = 0;
    *lptr = ptr+1;
    return atoi(line);
}

int
readMonth(char** lptr) {
    char* ptr = *lptr;
    char* line = *lptr;
    while (*ptr != ' ') ptr++;
    *ptr = 0;
    *lptr = ptr+1;
    int i;
    int m;
    for (i=0, m=1; i<12; i++, m++)
        if (strncmp(line, MONTHS[i], 3) == 0)
            return m;
            
    handle_error("invalid month");
    return -1;
}

int
countHits(const char* line, const char* key) {
    int n = strlen(key);
    int k = strlen(line) - n;
    int hits = 0;
    for (int i = 0; i < k; i++, line++) {
        if (*line == *key && memcmp(line, key, n) == 0) {
            ++hits;
            i += n;
        }
    }
    
    return hits;
}

void
writeTweet(char *tweet, const int fn, const int ln, const int hits,
           const int month, const int day, char* line) {
    short* ptr1 = (short*) tweet;
    *ptr1 = (short) fn;
    int* ptr2 = (int*) (tweet + 2);
    *ptr2 = ln;
    *(tweet+6) = (char) hits;
    *(tweet+7) = (char) month;
    *(tweet+8) = (char) day;
    
    int n = TSIZE-9;
    for (int i=strlen(line); i<n; i++) line[i] = ' '; // padding
    memcpy(tweet+9, line, n);
}

char
isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}


/**
 * MAIN PARALLEL EXECUTION
 */
void
parallel(const char* key, int rank, int size, int readedLines) {

    // Step 0 or 1
    // O(n-1)
    for(int op = 0, step = 0; op < size; op++, step ^= 1) {
        // Calculate who is sender and receiver
        int sender =    ((step == 0) && ( rank == 0 || rank % 2 == 0))
                        || (step == 1 && (rank != 0 && rank % 2 != 0 && !isLastProc(rank, size)));
                        
        int receiver = (step == 0 &&  rank != 0 && rank % 2 != 0)
                       || (step == 1 && ( rank != 0 && rank % 2 == 0 && !isLastProc(rank, size)));
                       
        ///////////////////////////////////////////////////
        // Compare HEIGH
        ///////////////////////////////////////////////////
        if(sender) {
            // 1-Send: worst Tweet to neighbors --------------------------------
            MPI_Send(
                TWEETS[step == 0 ? readedLines - 1 : 0], // Last/Worst Tweet
                TSIZE,
                MPI_CHAR,
                rank + 1,
                WORSTTWEET,
                MPI_COMM_WORLD);
#ifdef DEBUG
            printf("A: RANK %d SEND WORST TWEET TO %d\n", rank, rank + 1);
#endif
            // See A-Recv ------------------------------------------------------
            
            // A-Recv: NUM, BEST, BEST OF --------------------------------------
            char numBest[S_INT + TSIZE + TSIZE];
            MPI_Recv(
                numBest,
                S_INT + TSIZE + TSIZE,
                MPI_CHAR,
                rank + 1,
                NUMBEST,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
                
            int *better = (int*) numBest;
#ifdef DEBUG
            printf("B-Recv: RANK %d GET FROM %d : (%d)\n", rank, rank + 1, *better);
#endif
            
            if(*better > 0) {
                // First is already send, and if better < 0 so it's already okay
                int k = readedLines - *better -1;
                // TODO?
                k = k < 0 ? 0 : k;
                
                int sendTweetsUp = 0;
                int sendTweetsDown = 0;
                
                
                for(int i = readedLines - 1; k < readedLines -1; k++, i--)
                    if(compare(numBest + S_INT, TWEETS[step == 0 ? k: i]) == -1) ++sendTweetsUp;
                    else if(*better > 1 && compare(numBest + S_INT + TSIZE, TWEETS[step == 0 ? i: k]) == -1)
                        ++sendTweetsDown;
                    else continue;
                    
                    
                int sendTweets = (sendTweetsUp + sendTweetsDown + 1) / 2; // Round up
                
                if(sendTweets > 0) {

                    MPI_Sendrecv_replace(
                        TWEETS[step == 0 ? readedLines - sendTweets : 0],
                        sendTweets * TSIZE,
                        MPI_CHAR,
                        rank + 1,
                        2,
                        rank + 1,
                        3,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
#ifdef DEBUG
                    printf("C-D. RANK %d SEND RECV %d (%d)\n", rank, rank + 1, sendTweets);
#endif
                    
                    brutto_start = MPI_Wtime();
                    bitonicMerge(0, readedLines, step == 0 ? ASCENDING : DESCENDING);
                    brutto_end += MPI_Wtime() - brutto_start;
                } else {
                    MPI_Send(
                        NULL, // Last/Worst Tweet
                        0,
                        MPI_CHAR,
                        rank + 1,
                        2,
                        MPI_COMM_WORLD);
                }
            }
        }
        ///////////////////////////////////////////////////
        // Compare Low
        ///////////////////////////////////////////////////
        if(receiver) {
            // A-Recv: Worst Tweet----------------------------------------------
            char getData[TSIZE];
            MPI_Recv(
                getData,
                TSIZE,
                MPI_CHAR,
                rank - 1,
                WORSTTWEET,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
                
#ifdef DEBUG
            printf("A. RANK %d GET WORST TWEET FROM %d\n", rank, rank - 1);
#endif
            
            char numBest[S_INT + TSIZE + TSIZE];
            int *better = (int*) numBest;
            
            // Search better Tweets than the worst from left neighbor
            int betterTweets;
            if(step == 0) {
                betterTweets = readedLines - 1;
                for(int c; betterTweets > -1; betterTweets--) {
                    c = compare((void *)TWEETS[betterTweets], (void*) getData);
                    if(c == 1 || c == 0) break;
                }
                *better = readedLines -1 - betterTweets;
            } else {
                betterTweets = 0;
                for(int c; betterTweets < readedLines; betterTweets++) {
                    c = compare((void *)TWEETS[betterTweets], (void*) getData);
                    if(c == 1 || c == 0) break;
                }
                *better = betterTweets;
            }
            
            
            if(*better > 0) {
                memcpy(numBest + S_INT, TWEETS[step == 0 ? readedLines - 1 : 0], TSIZE); // Best
                if(*better > 1) {
                    memcpy(numBest + S_INT + TSIZE, TWEETS[step == 0 ? betterTweets + 1 : *better - 1], TSIZE); // Worst OF
                }
            }
            
            // SEND WORST AND BEST TWEET
            MPI_Send(
                numBest,
                S_INT + TSIZE + TSIZE,
                MPI_CHAR,
                rank - 1,
                NUMBEST,
                MPI_COMM_WORLD);
                
#ifdef DEBUG
            printf("B-Send: RANK %d SEND TO %d: (%d) \n", rank, rank - 1, *better);
#endif
            
            //------------------------------------------------------------------
            if(*better > 0) {
                MPI_Status status;
                int  nbytes = 0;
                
                MPI_Probe(
                    rank - 1,
                    2,
                    MPI_COMM_WORLD,
                    &status);

                if(MPI_Get_count(&status, MPI_CHAR, &nbytes) == MPI_UNDEFINED)
                    handle_error("ERROR: MPI_Get_count TAG 2");
                    
                int tweets = nbytes / TSIZE;

                if(tweets > 0) {
                    // Wir bekommenn die anderen Tweets zurueck
                    char *worstTweets = malloc(tweets * TSIZE);
                    if(worstTweets == NULL) handle_error("NO MEMORY FOR WORST TWEETS");
                    
                    MPI_Recv(
                        worstTweets,
                        tweets * TSIZE,
                        MPI_CHAR,
                        rank - 1,
                        2,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE);
#ifdef DEBUG
                    printf("C. RANK %d GET WORST %d TWEETS FROM %d\n",rank, tweets, rank  - 1);
#endif
                    MPI_Request r;
                    MPI_Send(
                        TWEETS[step == 0 ? readedLines - tweets - 1: 0],
                        tweets * TSIZE,
                        MPI_CHAR,
                        rank - 1,
                        3,
                        MPI_COMM_WORLD
                    );
#ifdef DEBUG
                    printf("D. RANK %d SEND BEST %d TWEETS TO %d\n", rank, tweets, rank  - 1);
#endif
                    memcpy(TWEETS[step == 0 ? readedLines - tweets: 0], worstTweets, tweets * TSIZE);
                    free(worstTweets);
                    
                    brutto_start = MPI_Wtime();
                    bitonicMerge(0, readedLines, step == 0 ? DESCENDING : ASCENDING);
                    brutto_end += MPI_Wtime() - brutto_start;
                }
            }
            
        }
    }
}

/**
 * MAIN EXECUTURE
 */
void
exec(const int numFiles, const int rank, const int size, const char* key) {

    FILE *file;
    
    char fileName[20];
    
    int startFile;
    if(numFiles > size)
        startFile = size / numFiles * rank; // 0
    else
        startFile = ((numFiles * 10 / size) * rank) / 10;
    
    int readFiles = (numFiles + size - 1) / size; // Roundup 2
    
    // Global Lines To Read
    linesToRead = TNUM * numFiles / size; // 48
    
    // Allocate size for tweets
    char *TWEETSDATA[readFiles];
    TWEETS = malloc(linesToRead * sizeof(char*)); // Array of pointers
    
    netto_start = MPI_Wtime();
    
    int iLine = 0;
    double startReadFile = MPI_Wtime();
    
    
    for(int f = startFile, fi = 0; f < readFiles + startFile; f++, fi++) {
    
        sprintf(fileName, FIN"%d",f);
        file = fopen(fileName, "r");
        if (file == NULL) handle_error("Cannon open File .\n");
        
        TWEETSDATA[fi] = malloc(linesToRead * TSIZE);
        char *data = TWEETSDATA[fi];
        if(data == NULL) handle_error("No more memory for TWEETSDATA.\n");
        
        int globalstart = 0;
        int globalend = 0;
        
        int myrank = -1;
        if(numFiles >= size){
           globalstart = 0;
           globalend = TNUM;
        } else {
            if(myrank != f)
            {
                globalstart = 0;
                
            }
            else{
                myrank = f;
                globalstart = rank * (TNUM / numFiles);
            }
            //globalstart = startFile * TNUM / (size / numFiles) ;
            
            //globalend   = globalstart + TNUM / (size / numFiles) -1;
            globalend   = globalstart + TNUM / (size / numFiles) -1;
        }
        // calc lines to read for last proc
        //if (isLastProc(rank, size)) {
        //    globalend = TNUM;
        //    linesToRead = TNUM - linesToRead * (size -1);
        //};
        
        printf("RANK %d file %d start %d end %d \n", rank, f,  globalstart, globalend);
        
        char buf[MAX_LINE_SIZE];
        char* line;
        
        int fn, ln, month, day, hits;
        
        for(int lines = 0, localLine = 0;
                lines < globalend && (line = fgets(buf, MAX_LINE_SIZE, file)) != NULL;
                ++lines) {
                
            if(lines < globalstart) continue;
            if(lines > globalend - 1) break;
            
            TWEETS[iLine] = data + (localLine * TSIZE);
            
            fn = readNumber(&line);
            
            ln = readNumber(&line);
            
            month = readMonth(&line);
            day = readNumber(&line);
            
            hits = countHits(line, key);
            
            writeTweet(TWEETS[iLine], fn, ln, hits, month, day, line);
            
            iLine++; // Global line
            localLine++;
        }
        
        // close file
        fclose(file);
        printf("readFiles %d\n", f);
    }
    
    brutto_start = MPI_Wtime();
    qsort(TWEETS[0], linesToRead, TSIZE, rank == 0 || rank % 2 == 0 ? compare : compareD);
    brutto_end = MPI_Wtime() - brutto_start;
    
    if(size > 1) parallel(key, rank, size, iLine);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    netto_end = MPI_Wtime() - netto_start;
    
    // Write tweets to file
    writeOrderedTweets(rank, TWEETS, linesToRead);
    
    // Free Space
    for(int i = 0; i < readFiles; i++) free(TWEETSDATA[i]);
    
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
    
    int numFiles = atoi(argv[2]);
    
    // Initialize MPI
    MPI_Init(NULL, NULL);
    
    // Get the number of processes and rank
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0) {
        printf("Number of process: '%d'\nFiles: '%d'\nKey: '%s'\n=== Start execution\n", size, numFiles, argv[1]);
    }
    
    // Execute main programm
    exec(numFiles, rank, size, argv[1]);

    if(rank == 0) {
        printf("\n Zeitmessung:  \n");
        printf("Rank: %d  - Netto ='%f'; Brutto ='%f'\n",rank, netto_end, brutto_end);
        
        double getTime[size][2];
        double send;
        
        MPI_Gather(&send, 1, MPI_DOUBLE, getTime, 2, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
        
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
