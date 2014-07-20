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

// #define TNUM 2400000 // Zeilen
#define FNUM 1 // Anzahl der Dateien
#define TSIZE 23
#define TNUM 2400000 // Zeilen

#define FIN "twitter.data."
//#define FIN "/mpidata/parsys14/gross/twitter.data.1"
#define U_MAX_BYTES 4

char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
int linesToRead;

double brutto_start, brutto_end, netto_start, netto_end;
///////////////////////////////////////////////////
// TWEETER STRUCT
///////////////////////////////////////////////////
typedef struct tweets {
    // Offset in file
    int offset;
    // File number
    unsigned int fn;
    // Line number
    unsigned int ln;
    // Number of hits
    unsigned int hits;
    
    char month;
    char day;
    char tweet[TSIZE];
} TDATA;

// Global Tweets array
TDATA **TWEETS;

/**
 * UTF-8. Return length of the UTF-8 Character.
// mask values for bit pattern of first byte in multi-byte
// UTF-8 sequences:
// 192 - 110xxxxx - for U+0080 to U+07FF
// 224 - 1110xxxx - for U+0800 to U+FFFF
// 240 - 11110xxx - for U+010000 to U+1FFFFF
static unsigned short mask[] = {192, 224, 240};

int u_getc(char chunk[], unsigned int offset, char bytes[]) {
    int i = 1;
    if ((chunk[offset] & mask[0]) == mask[0]) i++; // Ist Char ein UTF-8 zeichen?
    if ((chunk[offset] & mask[1]) == mask[1]) i++;
    if ((chunk[offset] & mask[2]) == mask[2]) i++;
    for(int j = 0; j < i; bytes[j] = chunk[offset + j++]); // Zeichen einlesen
    bytes[i] = '\0'; // Letztes Byte null setzten
    return i;
}
*/

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
void writeOrderedTweets(int rank, TDATA **T, int size) {
    char fileName[20];
    sprintf(fileName, "twitter.sort.%d", rank);
    
    remove(fileName);
    FILE* fu = fopen(fileName, "a");
    TDATA *t;
    for (int i = 0; i < size; i++) {
        t = T[i];
        fprintf(fu,  "R:%d FN:%d LN:%d - Hits: %d\n",rank,  t->fn, t->ln, t->hits);
    }
    
    fclose(fu);
}

/**
 * Swap Tweet
 */
void swap(int i, int j) {
    TDATA *t = TWEETS[i];
    TWEETS[i] = TWEETS[j];
    TWEETS[j] = t;
}

/**
 * Compare Tweet
 */
int compare(const void *ptr1, const void *ptr2) {

    TDATA *t1 = (TDATA*) ptr1;
    TDATA *t2 = (TDATA*) ptr2;
    
    if(t1->hits > t2->hits) return -1;
    if(t2->hits > t1->hits) return 1;
    
    if(t1->month > t2->month) return -1;
    if(t2->month > t1->month) return 1;
    
    if(t1->day > t2->day) return -1;
    if(t2->day > t1->day) return 1;
    
    for (int i = 0; i < TSIZE; i++) {
        if (t1->tweet[i] > t2->tweet[i]) return -1;
        if (t2->tweet[i] > t1->tweet[i]) return  1;
    }
    
    return 0;
}

/**
 * Bitonic sort of non power of two
 */
void bitonicMerge(int lo, int n, int dir) {
    if ( n > 1) {
        int m = n / 2;
        for (int i=lo; i<lo+m; i++) if(dir == compare(TWEETS[i], TWEETS[i + m])) swap(i, i + m);
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

/**
 * WRITE TWEET TO MEMORY
 */
void writeTweet(TDATA *tweet,
                const unsigned int fn, // File number
                const unsigned int ln, // Line numver
                const unsigned int hits, // Hits
                const unsigned int month, // Month
                const unsigned int day, // Day
                char* line, // Char line
                long offset // Offset in file
               ) {
    tweet->offset = offset;
    tweet->fn = fn;
    tweet->ln = ln;
    
    tweet->hits = hits;
    tweet->month = (char) month;
    tweet->day = (char) day;
    
    memset(tweet->tweet, ' ', TSIZE);
    int l = strlen(line);
    memcpy(tweet->tweet, line, l < TSIZE ? l: TSIZE);
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

/**
 * ALLOCALTE ONE BLOG MEMORY FOR SPECIFIED LINES
 */
TDATA **allocTweets(int lines) {
    TDATA **tweets = calloc(lines, sizeof(TDATA*)); // Array of pointers
    TDATA *data    = calloc(lines, sizeof(TDATA)); // The data
    
    for (int i = 0; i < lines ; i++) tweets[i] = &(data[i]);
    
    return tweets;
}

/**
 * GET TWEETS FROM FILE BY SPECIFIED FILE NUMBER AND OFFSET
 */
TDATA **getTweetFromFile(FILE *files[], const char* key,  unsigned int size, int data[][2]) {

    TDATA **tweetsFromFile = allocTweets(size);
    char* line;
    char buf[MAX_LINE_SIZE];
    unsigned int fn;
    unsigned int ln;
    unsigned int month;
    unsigned int day;
    unsigned int hits;
    
    for(int  i = 0; i < size; ++i) {
        // Seek to position
        fseek(files[data[i][0]], data[i][1], SEEK_SET);
        
        // Get line
        line = fgets(buf, MAX_LINE_SIZE, files[data[i][0]]);
        if(line == NULL) handle_error("file is null");
        
        fn = readNumber(&line);
        
        ln = readNumber(&line);
        
        month = readMonth(&line);
        day = readNumber(&line);
        
        hits = countHits(line, key);
        
        writeTweet(tweetsFromFile[i], data[i][0], ln, hits, month, day, line, data[i][1]);
    }
    
    return tweetsFromFile;
}

/**
 * COPY TWEET DATA FROM RIGHT TO LEFT
 */
void copyTweet(TDATA *t1, TDATA *t2) {
    t1->offset = t2->offset;
    t1->ln = t2->ln;
    t1->fn = t2->fn;
    t1->hits = t2->hits;
    t1->month = t2->month;
    t1->day = t2->day;
    memcpy(t1->tweet, t2->tweet, TSIZE);
}

/**
 * MAIN PARALLEL EXECUTION
 */
void parallel(FILE* files[], const char* key, int rank, int size, int readedLines) {
    if (rank == 0) {
        printf("Number of Processes spawned: %d\n", size);
    }
    
    // First we go down
    int down = 1;
    for(int next = size > 1; next > 0; ) {
    
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
            // Send last (worst) tweet to right
            printf("1.RANK %d SEND LAST TWEET TO %d\n",rank, rank + next);
            int sendData[2] = { TWEETS[linesToRead - 1]->fn, TWEETS[linesToRead - 1]->offset};
            MPI_Send(sendData, 2, MPI_INT, rank + next, 0, MPI_COMM_WORLD);
            
            MPI_Status status;
            int  nbytes = 0;
            MPI_Probe(rank + next, 2, MPI_COMM_WORLD, &status);
            if(MPI_Get_count(&status, MPI_INT, &nbytes) == MPI_UNDEFINED) {
                handle_error("Error msg");
            }
            
            // kopie von nbytes, da es nach MPI_Recv komische zahl hat
            int tweets = nbytes / 2;
            
            printf("5.RANK %d GET %d TWEET FROM %d\n",rank, tweets, rank + next);
            
            int getData[tweets][2];
            MPI_Recv(getData, tweets * 2, MPI_INT, rank + next, 2, MPI_COMM_WORLD, &status);
            
            if(tweets != 0) {
                // Get Tweets from file
                TDATA **tweetsFromRight = getTweetFromFile(files, key, tweets, getData);
                
                int j = linesToRead - 1; // Iterator goes from bad to good tweets
                int g = 0;
                
                // Welche Tweets sind  besser?
                int i = 0;
                for(; j > -1 && g < tweets; j--, i++) {
                    if(compare(tweetsFromRight[g], TWEETS[j]) == -1) {
                        g++;
                    } else {
                        break;
                    }
                }
                
                int h = linesToRead - i; // Iterator goes from bad tweets - tweets that are better
                int sendData[i][2];
                
                for(int k = 0; k  < i ; k++, h++) {
                    sendData[k][0] = TWEETS[h]->fn;
                    sendData[k][1] = TWEETS[h]->offset;
                    copyTweet(TWEETS[h], tweetsFromRight[k]);
                }
                
                free(tweetsFromRight[0]);
                free(tweetsFromRight);
                
                printf("6.RANK %d SEND %d TWEET TO %d\n",rank, i, rank + next);
                
                // Unblocked!
                MPI_Request request;
                MPI_Send(&sendData, i * 2, MPI_INT, rank + next, 3, MPI_COMM_WORLD);
            }
        }
        ///////////////////////////////////////////////////
        // Compare Low
        ///////////////////////////////////////////////////
        if(receiver) {
            int getData[1][2];
            printf("2.RANK %d GET LAST TWEET FROM %d\n", rank, rank - next);
            MPI_Recv(getData, 2, MPI_INT, rank - next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            TDATA **tweetsFromFile = getTweetFromFile(files, key, 1, getData);
            int betterTweets = 0;
            for(int i = 0; betterTweets < linesToRead; i++) {
                int c = compare( TWEETS[i], tweetsFromFile[0] );
                if( c == -1 )
                    betterTweets++;
                else
                    break;
            }
            
            if(betterTweets > 0) {
                int sendData[betterTweets][2];
                
                for(int i = 0; i < betterTweets; i++) {
                    sendData[i][0] = TWEETS[i]->fn;
                    sendData[i][1] = TWEETS[i]->offset;
                };
                
                MPI_Send(sendData, betterTweets * 2, MPI_INT, rank - next, 2, MPI_COMM_WORLD);
                printf("3.RANK %d SEND %d TWEET TO %d\n",rank, betterTweets, rank - next);
                
                MPI_Status status;
                int  nbytes = 0;
                MPI_Probe(rank - next, 3, MPI_COMM_WORLD, &status);
                if(MPI_Get_count(&status, MPI_INT, &nbytes) == MPI_UNDEFINED) {
                    handle_error("ERROR: MPI_Get_count TAG 3");
                }
                
                // kopie von nbytes, da es nach MPI_Recv komische zahl hat
                int tweets = nbytes  / 2;
                
                int getData[tweets][2];
                // Wir bekommenn die anderen Tweets zurueck
                MPI_Recv(getData, tweets * 2, MPI_INT, rank - next, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("4.RANK %d GET %d TWEET FROM %d\n",rank, tweets, rank  - next);
                
                TDATA **tweetsFromLeft = getTweetFromFile(files, key, tweets, getData);
                
                // writeOrderedTweets(tweetsFromLeft, j);
                for(int k = 0; k < tweets; k++) {
                    copyTweet(TWEETS[k], tweetsFromLeft[k]);
                }
                
                free(tweetsFromLeft[0]);
                free(tweetsFromLeft);
                
            } else {
                MPI_Send(NULL, 0, MPI_INT, rank - next, 2, MPI_COMM_WORLD);
            }
        }
        
        if(next > size / 2 ) {
            // Go UP!
            down = 0;
        }
        
        if(down) next++;
        else next--;
        
        bitonicSort(0, readedLines, ASCENDING);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
}

void exec(const int numFiles, const int rank, const int size, const char* key) {
    FILE *files[FNUM];
    
    char fileName[20];
    
    linesToRead = TNUM / size;
    //printf("Lines %d \n", linesToRead);
    // Allocate size
    TWEETS = allocTweets(FNUM * linesToRead * size);
    int iLine = 0;
    for(int f = 0; f < FNUM; f++) {
        sprintf(fileName, FIN"%d",f);
        files[f] = fopen(fileName, "rb");
        if (files[f] == NULL) handle_error("Cannon open File .\n");
        
        int globalstart = rank * linesToRead;
        int globalend   = globalstart + linesToRead;
        int sortStart = iLine;
        
        // Der letzter Prozessor
        if (isLastProc(rank, size)) {
            globalend = TNUM;
            linesToRead = TNUM - linesToRead * (size -1);
        };
        
        char buf[MAX_LINE_SIZE];
        
        unsigned int fn;
        unsigned int ln;
        unsigned int month;
        unsigned int day;
        unsigned int hits;
        int offset;
        char* line;
        
        // printf("iLines %d \n", iLine);
        for(int lines = 0; iLine < linesToRead * FNUM && (line = fgets(buf, MAX_LINE_SIZE, files[f])) != NULL; ++lines) {
            if(lines < globalstart) continue;
            if(lines > globalend -1) break;
            
            offset = ftell(files[f]) - strlen(line);
            
            fn = readNumber(&line);
            
            ln = readNumber(&line);
            
            month = readMonth(&line);
            day = readNumber(&line);
            
            hits = countHits(line, key);
            
            writeTweet(TWEETS[iLine], f, ln, hits, month, day, line, offset);
            
            iLine++;
        }
        
        
        // 1. Schritt Sortiere Locale Tweets
        //bitonic(iLine);
        //qsort(TWEETS, iLine, sizeof(TDATA*), compare);
        brutto_start = MPI_Wtime();
        bitonicSort(sortStart, iLine, ASCENDING);
        
        // Warten bis alle Prozessoren hier sind
        MPI_Barrier(MPI_COMM_WORLD);
        //writeOrderedTweets(rank, TWEETS, linesToRead, 1);
        
        // 2. Schritt: Tausche Tweets aus
        parallel(files, key, rank, size, iLine);
        brutto_end = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    //qsort(TWEETS, linesToRead, sizeof(TDATA*), compare);
    bitonicSort(0, linesToRead, ASCENDING);
    
    writeOrderedTweets(rank, TWEETS, linesToRead);
    
    free(TWEETS);
    
    for(int i = 0; i < FNUM; i++) fclose(files[i]);
}



int main(int argc, char** argv) {
    
    if(argc != 2) {
        fprintf(stderr, "Please specify search key");
        return EXIT_FAILURE;
    }
    
    // Initialize MPI
    MPI_Init(NULL, NULL);
    netto_start = MPI_Wtime();
    // Get the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Execute main programm
    exec(FNUM, rank, size, argv[1]);
    
    
    char processor_name[MPI_MAX_PROCESSOR_NAME]; // Get the name of the processor
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Print off a hello world message
    // printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    netto_end = MPI_Wtime();
    printf("Netto: %f seconds from Rank %d\n", netto_end-netto_start, rank);
    printf("Brutto: %f seconds from Rank %d\n", brutto_end-brutto_start, rank);
    MPI_Finalize();
    return EXIT_SUCCESS;
}


