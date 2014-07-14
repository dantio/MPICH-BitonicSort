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
#define TSIZE 32
#define TNUM 8 // Zeilen
#define FIN "twitter.data10"
#define U_MAX_BYTES 4

char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

int linesToRead;
char **TWEETS; // Pointers to Array of pointers

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

void printTweet(const char* t, char* buffer) {
    int i;
    int n = 1;
    buffer[0] = '[';
    for (i=0; i < TSIZE; i++) {
        int k = t[i];
        int l = sprintf(buffer + n, "%2x ", k<0 ? k+256 :k);
        n +=l;
    }
    buffer[n] = ']';
    buffer[n+1] = '\n';
    buffer[n+2] = '\0';
    printf("%s", buffer);
}

void writeOrderedTweets() {
    char buffer[linesToRead * TSIZE * 3];
    for (int i = 0; i < linesToRead; i++) {
        printTweet(TWEETS[i], buffer + i * TSIZE);
    }
    
}

void swap(int i, int j) {
    char *t = TWEETS[i];
    TWEETS[i] = TWEETS[j];
    TWEETS[j] = t;
}


int compare(const char t1[TSIZE], const char t2[TSIZE]) {
    for (int i = 6; i < TSIZE; i++) {
        if (t1[i] > t2[i]) return -1;
        if (t2[i] > t1[i]) return  1;
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


int readNumber(char** lptr) {
    char* ptr = *lptr;
    char* line = *lptr;
    while (*ptr != ' ') ptr++;
    *ptr = 0;
    *lptr = ptr+1;
    return atoi(line);
}

int readMonth(char** lptr) {
    char* ptr = *lptr;
    char* line = *lptr;
    while (*ptr != ' ') ptr++;
    *ptr = 0;
    *lptr = ptr+1;
    int i, m;
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


void writeTweet(int i, const int fn, const int ln, const int hits,
                const int month, const int day, char* line) {
    TWEETS[i][0] = fn;
    TWEETS[i][2] = ln;
    
    TWEETS[i][6] = (char) hits;
    TWEETS[i][7] = (char) month;
    TWEETS[i][8] = (char) day;
    int n = TSIZE - 9;
    char *ptr = &TWEETS[i][9];
    for (int j = strlen(line); i < n; i++) line[j] = ' '; // padding
    memcpy(ptr, line, n);
    
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

char **allocTweets(int lines) {
    char **tweets = ( char** ) malloc(lines * sizeof( char* ));
    for (int i = 0; i < lines; i++ ) {
        if (( tweets[i] = ( char* ) malloc( TSIZE )) == NULL)
            handle_error("NO SPACE");
    }
    return tweets;
}

void exec(const int rank, const int size, const char* key) {

    FILE* file = fopen(FIN, "r");
    if (file == NULL) handle_error("Cannon open File .\n");
    
    int BUFFER_SIZE = MAX_LINE_SIZE;
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
    
    char buf[BUFFER_SIZE];
    
    unsigned int lines = 0, fn, ln, month, day, hits;
    char* line;
    for(int i = 0; (line = fgets(buf, BUFFER_SIZE, file)) != NULL; ++lines) {
        if(lines < globalstart) continue;
        if(lines > globalend -1) break;
        
        fn = readNumber(&line);
        
        ln = readNumber(&line);
        
        month = readMonth(&line);
        day = readNumber(&line);
        
        hits = countHits(line, key);
        
        writeTweet(i, fn, ln, hits, month, day, line);
        
        i++;
    }
    
    fclose(file);
    
    bitonic(linesToRead);
    
    // Warten bis alle Prozessoren hier sind
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(!isLastProc(rank, size) && (rank == 0 ||rank % 2 == 0)) {
        // Sende letzten TWEET
        MPI_Send(TWEETS[linesToRead -1], TSIZE, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
    }
    
    if(rank != 0 && rank % 2 != 0) {
        char lastTweet[TSIZE];
        MPI_Recv(&lastTweet, TSIZE, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Its der letzte TWEET groeser als die ersten im rechten Knoten?
        int j = 0;
        for(; j < linesToRead; j++) {
            if(compare(&lastTweet, TWEETS[j]) == 1) {
                continue;
            } else break;
        }
        
        // nbytes ist hier 0 rechter Knoten keine kleineren hat es nicht so ist
        int nbytes = j * TSIZE;
        MPI_Send(TWEETS[0], nbytes, MPI_CHAR, rank - 1, 2, MPI_COMM_WORLD);
        
        // Wir bekommenn die anderen Tweets zurueck
        if(nbytes != 0) {
            MPI_Recv(TWEETS[0], nbytes, MPI_CHAR, rank - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    
    if(!isLastProc(rank, size) && (rank == 0 ||rank % 2 == 0)) {
    
        MPI_Status status;
        int  nbytes = 0;
        // Wieviel will der rechte knoten uns schicken?
        MPI_Probe(rank + 1, 2, MPI_COMM_WORLD, &status);
        // It will block the caller until a message is ready
        
        if(MPI_Get_count(&status, MPI_CHAR, &nbytes) == MPI_UNDEFINED) {
            handle_error("Error msg");
        }
        
        // kopie von nbytes, da es nach MPI_Recv komische zahl hat
        int bytes = nbytes;
        int lines = bytes / TSIZE;
        
        // Platz fuer die Tweets aus dem rechten Knoten
        char **movedTweets = allocTweets(lines);
        
        // Finally, receive the message with a correctly sized buffer...
        MPI_Recv(movedTweets[0], nbytes, MPI_CHAR, rank + 1, 2, MPI_COMM_WORLD, &status);
        if(nbytes != 0) {
        
        
            MPI_Send(TWEETS[lines - 1], bytes, MPI_CHAR, rank + 1, 3, MPI_COMM_WORLD);
            for(int k = 0; lines - 1 < linesToRead; lines++, k++) {
                char *freeMee = TWEETS[lines];
                TWEETS[lines - 1] = movedTweets[k];
                free(freeMee);
            }
        }
    }
    
    bitonic(linesToRead);
    MPI_Barrier(MPI_COMM_WORLD);
    writeOrderedTweets();
    //printf("---------------\n");
    
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
    
    char processor_name[MPI_MAX_PROCESSOR_NAME]; // Get the name of the processor
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Print off a hello world message
    // printf("Hello world from processor %s, rank %d out of %d processors\n",
    //      processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return EXIT_SUCCESS;
}


