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
#define MAX_LINE_SIZE 1000

// #define TNUM 2400000 // Zeilen
#define FNUM 16	
#define TSIZE 24
#define TNUM 512 // Zeilen

#define FIN "twitter.data1000"
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


void writeOrderedTweets(int rank, TDATA **T, int size) {
    char buffer[500] = "";
    char allTweets[5000] = "";
    sprintf(buffer, "-----  Rank %d ------ \n", rank);
    strcat(allTweets, buffer);
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
       return 0;
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
    
    for (int k = strlen(line); k < TSIZE; k++) line[k] = ' '; // padding
    memcpy(tweet->tweet, line, TSIZE);
    
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

TDATA **allocTweets(int lines) {
    TDATA *data = ( TDATA* ) malloc(lines * sizeof(TDATA));
    TDATA **tweets = ( TDATA** ) malloc(lines * sizeof(TDATA*));
    for (int i = 0; i < lines; i++) tweets[i] = &(data[i]);
    
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

void parallel(FILE* file, const char* key, int rank, int size) {
    if(!isLastProc(rank, size) && (rank == 0 ||rank % 2 == 0)) {
        // Sende letzten
        MPI_Send(&(TWEETS[linesToRead - 1]->offset), 1, MPI_LONG, rank + 1, 0, MPI_COMM_WORLD);
        printf("1.RANK %d SEND LAST TWEET\n",rank);
    }
    
    if(rank != 0 && rank % 2 != 0) {
    
        long lastTweetLn[1];
        MPI_Recv(lastTweetLn, 1, MPI_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("2.RANK %d GET LAST TWEET\n",rank);
        
        // Ist der letzte TWEET groeser als die ersten im rechten Knoten?
        int j = 0;
        
        TDATA **tweetsFromFile = getTweetFromFile(file, key, 1, lastTweetLn);
        //writeOrderedTweets(tweetsFromFile, 1);
        
        for(int i = 0; j < linesToRead; i++) {
            int c = compare( TWEETS[i],tweetsFromFile[0] );
            if( c == -1 )
                j++;
            else
                break;
        }
        
        
        if(j > 0) {
            long sendTweetsLn[j];
            
            for(int i = 0; i < j; i++) sendTweetsLn[i] = TWEETS[i]->offset;
            // printf("-- %d TWEET\n ", TWEETS[2]->ln);
            
            MPI_Send(sendTweetsLn, j, MPI_LONG, rank - 1, 2, MPI_COMM_WORLD);
            printf("3.RANK %d SEND %d TWEET\n",rank, j);
            
            // Wieviel will der linke knoten uns schicken?
            MPI_Status status;
            int  nbytes = 0;
            // Wieviel will der rechte knoten uns schicken?
            MPI_Probe(rank - 1, 3, MPI_COMM_WORLD, &status);
            // It will block the caller until a message is ready
            
            if(MPI_Get_count(&status, MPI_LONG, &nbytes) == MPI_UNDEFINED) {
                handle_error("Error msg");
            }
            
            // kopie von nbytes, da es nach MPI_Recv komische zahl hat
            int tweets = nbytes;
            
            long recTweets[tweets];
            // Wir bekommenn die anderen Tweets zurueck
            MPI_Recv(recTweets, tweets, MPI_LONG, rank - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("4.RANK %d GET %d TWEET\n",rank, tweets);
            
            TDATA **tweetsFromLeft = getTweetFromFile(file, key, tweets, recTweets);
            
            // writeOrderedTweets(tweetsFromLeft, j);
            for(int k = 0; k < tweets; k++) {
                //TDATA *freeMee = TWEETS[k];
                TWEETS[k] = tweetsFromLeft[k];
                //free(freeMee);
            }
            
            
        } else {
            MPI_Send(NULL, 0, MPI_LONG, rank - 1, 2, MPI_COMM_WORLD);
        }
    }
    
    
    if(!isLastProc(rank, size) && (rank == 0 || rank % 2 == 0)) {
    
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
        
        printf("5.RANK %d GET %d TWEET\n",rank, tweets);
        
        long sendTweetsLn[tweets];
        // Finally, receive the message with a correctly sized buffer...
        MPI_Recv(&sendTweetsLn, tweets, MPI_LONG, rank + 1, 2, MPI_COMM_WORLD, &status);
        
        if(tweets != 0) {
        
            // Platz fuer die Tweets aus dem rechten Knoten
            TDATA **tweetsFromRight = getTweetFromFile(file, key, tweets, sendTweetsLn);
            int j = 0;
            int k = 0;
            
            int h = 0;
            // Welche Tweets sind  besser?
            int i = 0;
            for(j = linesToRead - 1; j > -1 && i < tweets; j--, i++) {
                if(compare(tweetsFromRight[tweets - 1], TWEETS[j]) == -1){
                //printf("xxxxxxxxx\n");
                continue;
                }
                else break;
            }
            
            h = linesToRead - i;
            printf("JJJ %d\n", j);
            
            
            for(; k  < i ; k++, h++) {
                //printf("%d - %d \n",TWEETS[linesToRead - tweets + j]->hits, tweetsFromRight[k]->hits);
                sendTweetsLn[k] = TWEETS[h]->offset;
                //TDATA *freeMee = TWEETS[linesToRead - tweets + j];
                TWEETS[h] = tweetsFromRight[k];
                //free(freeMee);
                
            }
            
            printf("6.RANK %d SEND %d TWEET\n",rank, k);
            
            // Unblocked!
            MPI_Send(&sendTweetsLn, i, MPI_LONG, rank + 1, 3, MPI_COMM_WORLD);
        }
    }
    
    
}

void exec(const int rank, const int size, const char* key) {

    FILE* file = fopen(FIN, "rb");
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
    writeOrderedTweets(rank, TWEETS, linesToRead);
    
    //parallel(file, key, rank, size);
    
    bitonic(linesToRead);
    
    parallel(file, key, rank, size);
    
    bitonic(linesToRead);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    writeOrderedTweets(rank, TWEETS, linesToRead);
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
    
    /*TDATA **tweets = allocTweets(300); 24476
    char line[32] = "hello";
    for(int i = 0; i < 300; i++)
      writeTweet(tweets, i, 2, 3, 4, 5, 6, line, 7);*/
    
    char processor_name[MPI_MAX_PROCESSOR_NAME]; // Get the name of the processor
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Print off a hello world message
    //  printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return EXIT_SUCCESS;
}


