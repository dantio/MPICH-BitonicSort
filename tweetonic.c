/*
 ============================================================================
 Name        : Tweetonic.c
 Author      : Daniil Tomilow
 ============================================================================
 */

#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


// 1  +  1  +  N   + 1  +   3  +    1  +   2  +  1 +  140 * 3 + 1
// 0  + ' ' + '10' + ' ' + 'MAR' + ' ' + '02' + ' '+ 'tweet' + '\n'
#define MAX_LINE_SIZE 722

// #define TNUM 2400000 // Zeilen
#define TSIZE 24
#define TNUM 8 // Zeilen

#define FIN "twitter.data10"
#define U_MAX_BYTES 4

char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

int linesToRead;

// = 66b
typedef struct tweets {
    long offset; // Offset in file 16b
    unsigned int fn; // File number 8b
    unsigned int ln; // Line number 8b
    unsigned int hits; // Hits 8b
    
    char month; // 1b
    char day; //1b
    char tweet[TSIZE]; // 24b
} TDATA;

TDATA **TWEETS; // Pointers to Array of pointers

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

char isPowerOfTwo(unsigned const int number) {
    return (
               (number != 0)
               &&
               ((number & (number - 1)) == 0)
           );
}


void writeOrderedTweets(TDATA **T, int size) {
    char buffer[500] = "";
    char allTweets[5000] = "";
    for (int i = 0; i < size; i++) {
        TDATA *t = T[i];
        sprintf(buffer, "%d - %d\n", t->ln, t->hits);
        strcat(allTweets, buffer);
    }
    //allTweets[strlen(allTweets)] = '\0';
    printf("%s \n", allTweets);
    
}

void swap(int i, int j) {
    TDATA *t = TWEETS[i];
    TWEETS[i] = TWEETS[j];
    TWEETS[j] = t;
}


int compare(const TDATA *t1, const TDATA *t2) {

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

void bitonic(unsigned int lines) {
    unsigned int i,j,ij,k,c;
    for (k = 2; k <= lines; k = 2 * k)
        for (j = k >> 1; j > 0; j = j >> 1)
            for (i = 0; i < lines; i++)
                if ((ij = i ^ j) > i && (c = compare(TWEETS[i], TWEETS[ij])))
                    if ((!(i & k) && c == 1) || ( (i & k) && c == -1))
                        swap(i,ij);
}


void handle_error(const char* msg) {
    MPI_Finalize();
    fprintf(stderr, "Error: %s\n", msg);
    exit(255);
}


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

int countHits(const char* line, const char* key) {
    int n = strlen(key);
    int k = strlen(line) - n;
    int i;
    int hits = 0;
    for (i = 0; i < k; i++, line++) {
        if (*line == *key) {
            if (memcmp(line, key, n) == 0)
                hits++;
        }
    }
    return hits;
}


void writeTweet(TDATA **T, unsigned int i, const unsigned int fn, const unsigned int ln, const unsigned int hits,
                const unsigned int month, const unsigned int day, char* line, long offset) {
    TDATA *tweet = T[i];
    
    tweet->offset = offset;
    tweet->fn = fn;
    tweet->ln = ln;
    
    tweet->hits = hits;
    tweet->month = (char) month;
    tweet->day = (char) day;
    
    int j = strlen(line);
    for (int k = 0; k < TSIZE; k++) line[k] = ' '; // padding
    memcpy(tweet->tweet, line, TSIZE);
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

TDATA **allocTweets(int lines) {
    // TODO
    int size = sizeof(TDATA);
    TDATA *data = ( TDATA* ) malloc(lines * size);
    TDATA **tweets = ( TDATA** ) malloc(lines * sizeof(TDATA*));
    for (int i = 0, j = 0; j < lines; i+= size, j++) tweets[j] = data + i;
    
    return tweets;
}

TDATA **getTweetFromFile(FILE *file, const char* key,  unsigned int size, long offset[]) {
    // TODO
    // Allocate size
    TDATA **tweetsFromFile = allocTweets(size);
    char* line;
    char buf[MAX_LINE_SIZE];
    unsigned int fn;
    unsigned int ln;
    unsigned int month;
    unsigned int day;
    unsigned int hits;
    
    for(int  i = 0; i < size; ++i) {
        fseek(file, offset[i], SEEK_SET);
        
        line = fgets(buf, MAX_LINE_SIZE, file);
        if(line == NULL) handle_error("file is null");
        //printf("%s", line);
        
        fn = readNumber(&line);
        
        ln = readNumber(&line);
        
        month = readMonth(&line);
        day = readNumber(&line);
        
        hits = countHits(line, key);
        
        writeTweet(tweetsFromFile, i, fn, ln, hits, month, day, line, offset[i]);
    }
    
    
    return tweetsFromFile;
}

void parallel(FILE* file, const char* key, int rank, int size){
    if(!isLastProc(rank, size) && (rank == 0 ||rank % 2 == 0)) {
        // Sende letzten
        // printf(">>>>%d  - %d\n",TWEETS[linesToRead - 1]->ln,  TWEETS[linesToRead - 1]->hits);
        MPI_Send(&(TWEETS[linesToRead - 1]->offset), 1, MPI_LONG, rank + 1, 0, MPI_COMM_WORLD);
    }
    
    if(rank != 0 && rank % 2 != 0) {
    
        long lastTweetLn[1];
        MPI_Recv(lastTweetLn, 1, MPI_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Ist der letzte TWEET groeser als die ersten im rechten Knoten?
        int j = 0;
        
        TDATA **tweetsFromFile = getTweetFromFile(file, key, 1, lastTweetLn);
        //writeOrderedTweets(tweetsFromFile, 1);
        
        for(int i = 0; j < linesToRead; i++) {
            if(compare( tweetsFromFile[0], TWEETS[i]) == 1)
                j++;
            else
                break;
        }
        
        printf("SEND %d tweets\n ", j);
        
        if(j > 0) {
            long sendTweetsLn[j];
            
            for(int i = 0; i < j; i++) sendTweetsLn[i] = TWEETS[i]->offset;
           // printf("-- %d TWEET\n ", TWEETS[2]->ln);
            
            MPI_Send(sendTweetsLn, j, MPI_LONG, rank - 1, 2, MPI_COMM_WORLD);
            
            // Wir bekommenn die anderen Tweets zurueck
            MPI_Recv(sendTweetsLn, j, MPI_LONG, rank - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            TDATA **tweetsFromLeft = getTweetFromFile(file, key, j, sendTweetsLn);
            // writeOrderedTweets(tweetsFromLeft, j);
            for(int k = 0; k < j; k++) {
                //TDATA *freeMee = &TWEETS[k];
                TWEETS[k] = tweetsFromLeft[k];
                //free(freeMee);
            }
            
        } else {
            MPI_Send(NULL, 0, MPI_LONG, rank - 1, 2, MPI_COMM_WORLD);
        }
    }
    
    
    if(!isLastProc(rank, size) && (rank == 0 ||rank % 2 == 0)) {
    
        MPI_Status status;
        int  nbytes = 0;
        // Wieviel will der rechte knoten uns schicken?
        MPI_Probe(rank + 1, 2, MPI_COMM_WORLD, &status);
        // It will block the caller until a message is ready
        
        if(MPI_Get_count(&status, MPI_LONG, &nbytes) == MPI_UNDEFINED) {
            handle_error("Error msg");
        }
        
        // kopie von nbytes, da es nach MPI_Recv komische zahl hat
        int tweets = nbytes;
        
        printf("RANK %d GET %d TWEET\n",rank, tweets);
        
        long sendTweetsLn[tweets];
        // Finally, receive the message with a correctly sized buffer...
        MPI_Recv(&sendTweetsLn, tweets, MPI_LONG, rank + 1, 2, MPI_COMM_WORLD, &status);
        
        char buff[200];
        
        if(tweets != 0) {
            // Platz fuer die Tweets aus dem rechten Knoten
            TDATA **tweetsFromRight = getTweetFromFile(file, key, tweets, sendTweetsLn);
           // writeOrderedTweets(tweetsFromRight, tweets);
            for(int i = linesToRead - tweets, j = 0; i < linesToRead; i++, j++) sendTweetsLn[j] = TWEETS[i]->offset;
            
            MPI_Send(&sendTweetsLn, tweets, MPI_LONG, rank + 1, 3, MPI_COMM_WORLD);
            
            for(int k = 0, m = linesToRead - tweets; m < linesToRead; m++, k++) {
                //TDATA *freeMee = TWEETS[m];
                TWEETS[m] = tweetsFromRight[k];
                //free(freeMee);
            }
        }
    }


}

void exec(const int rank, const int size, const char* key) {

    FILE* file = fopen(FIN, "r");
    if (file == NULL) handle_error("Cannon open File .\n");
    
    linesToRead = TNUM / size;
    
    int globalstart = rank * linesToRead;
    int globalend   = globalstart + linesToRead;
    
    // Der letzter Prozessor
    if (isLastProc(rank, size)) {
        globalend = TNUM;
        linesToRead = TNUM - linesToRead * (size -1);
    };
    
    //printf("Lines: %d\n", linesToRead);
    
    // Allocate size
    TWEETS = allocTweets(linesToRead);
    
    char buf[MAX_LINE_SIZE];
    
    unsigned int fn;
    unsigned int ln;
    unsigned int month;
    unsigned int day;
    unsigned int hits;
    long offset;
    char* line;
    for(int i = 0, lines = 0; (line = fgets(buf, MAX_LINE_SIZE, file)) != NULL; ++lines) {
        if(lines < globalstart) continue;
        if(lines > globalend -1) break;
        
        offset = ftell(file) - strlen(line);
        
        fn = readNumber(&line);
        
        ln = readNumber(&line);
        
        month = readMonth(&line);
        day = readNumber(&line);
        
        hits = countHits(line, key);
        
        writeTweet(TWEETS, i, fn, ln, hits, month, day, line, offset);
        
        i++;
    }
    
    
    
    bitonic(linesToRead);
    
    // Warten bis alle Prozessoren hier sind
    MPI_Barrier(MPI_COMM_WORLD);
    writeOrderedTweets(TWEETS, linesToRead);
    
	parallel(file, key, rank, size);
    
    bitonic(linesToRead);
    
    parallel(file, key, rank, size);
    
    bitonic(linesToRead);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
     writeOrderedTweets(TWEETS, linesToRead);
    //printf("---------------\n");
    fclose(file);
    // FREE
}



int main(int argc, char** argv) {
    if(argc != 2) {
        fprintf(stderr, "Please specify search key");
        return EXIT_FAILURE;
    }
    
    // Initialize MPI
    MPI_Init(NULL, NULL);
    
    // Get the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    exec(rank, size, argv[1]);
    
    // writeTweet(t, 1, 2, 3, 4, 5, 6, line, 7);
    char processor_name[MPI_MAX_PROCESSOR_NAME]; // Get the name of the processor
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Print off a hello world message
    //  printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return EXIT_SUCCESS;
}


