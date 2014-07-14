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

typedef struct tweets{
	unsigned int fn;
	unsigned int ln;
	unsigned int hits;
	unsigned int month;
	unsigned int day;
	char tweet[TSIZE];
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


void writeOrderedTweets() {
    char buffer[500];
    char allTweets[5000];
    for (int i = 0; i < linesToRead; i++) {
        TDATA *t = TWEETS[i];
        sprintf(buffer, "%d - %d %d %d\n", i, t->fn, t->ln,t->hits);
        strcat(allTweets, buffer);
    }
    allTweets[strlen(allTweets)] = '\0';
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


void writeTweet(unsigned int i, const unsigned int fn, const unsigned int ln, const unsigned int hits,
                const unsigned int month, const unsigned int day, char* line) {
    TDATA *tweet = TWEETS[i];
    tweet->fn = fn;
    tweet->ln = ln;

    tweet->hits = hits;    
    tweet->month = month;
    tweet->day = day;
    
    int j = strlen(line);
    for (int k = 0; k < TSIZE; k++) line[k] = ' '; // padding
      memcpy(tweet->tweet, line, TSIZE);
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

TDATA **allocTweets(int lines) {
    TDATA *data = ( TDATA* ) malloc(lines * (sizeof(TDATA) + TSIZE));
    TDATA **tweets = ( TDATA** ) malloc(lines * sizeof(TDATA*));
    for (int i = 0; i < lines; ++i ) {
        tweets[i] = &(data[lines * i]);
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
    
    unsigned int fn;
    unsigned int ln;
    unsigned int month;
    unsigned int day;
    unsigned int hits;
    char* line;
    for(int i = 0, lines = 0; (line = fgets(buf, BUFFER_SIZE, file)) != NULL; ++lines) {
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
    MPI_Barrier(MPI_COMM_WORLD);/*
    
    if(!isLastProc(rank, size) && (rank == 0 ||rank % 2 == 0)) {
        // Sende letzten TWEET
        MPI_Send(&(TWEETS[linesToRead - 1][0]), TSIZE, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
    }
    
    if(rank != 0 && rank % 2 != 0) {
        char lastTweet[TSIZE];
        MPI_Recv(&lastTweet, TSIZE, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Ist der letzte TWEET groeser als die ersten im rechten Knoten?
        int j = 0;
        for(; j < linesToRead; j++) {
            if(compare(&lastTweet, TWEETS[j]) == 1) {
                continue;
            } else break;
        }
        
        if(j > 0) {
            
            // nbytes ist hier 0 rechter Knoten keine kleineren hat es nicht so ist
            int nbytes = j * TSIZE;
            char b = 'b';
            MPI_Send(&b, 1, MPI_CHAR, rank - 1, 2, MPI_COMM_WORLD);
            
            
            // Wir bekommenn die anderen Tweets zurueck
            if(nbytes != 0) {
                TDATA **movedTweets = allocTweets(j);
                
                MPI_Recv(&(movedTweets[0][0]), nbytes, MPI_CHAR, rank - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int k = 0; k < j; k++) {
                    TDATA *freeMee = TWEETS[k];
                    TWEETS[k] = movedTweets[k];
                    free(freeMee);
                }
            }
        } else {
            MPI_Send(NULL, 0, MPI_CHAR, rank - 1, 2, MPI_COMM_WORLD);
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
        TDATA **movedTweets = allocTweets(lines);
        printf("send %d\n", bytes);
        // Finally, receive the message with a correctly sized buffer...
        MPI_Recv(&(movedTweets[0][0]), bytes, MPI_CHAR, rank + 1, 2, MPI_COMM_WORLD, &status);
        
        char buff[200];
        
        printTweet(movedTweets[0], buff);
        
        if(bytes != 0) {
            MPI_Send(&(TWEETS[linesToRead - lines][0]), bytes, MPI_CHAR, rank + 1, 3, MPI_COMM_WORLD);
            for(int k = 0, m = linesToRead - lines; m < linesToRead; m++, k++) {
                TDATA *freeMee = TWEETS[m];
                TWEETS[m] = movedTweets[k];
                free(freeMee);
            }
        }
    }
    
    bitonic(linesToRead);
    MPI_Barrier(MPI_COMM_WORLD);*/
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


