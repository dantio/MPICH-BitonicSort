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
#define TNUM 2 // Zeilen
#define FIN "twitter.data10"

// mask values for bit pattern of first byte in multi-byte
// UTF-8 sequences:
// 192 - 110xxxxx - for U+0080 to U+07FF
// 224 - 1110xxxx - for U+0800 to U+FFFF
// 240 - 11110xxx - for U+010000 to U+1FFFFF
#define U_MAX_BYTES 4
char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
int linesToRead;
static unsigned short mask[] = {192, 224, 240};

/**
 * UTF-8. Return length of the UTF-8 Character.
 */
int u_getc(char chunk[], unsigned int offset, char bytes[]) {
    int i = 1;
    if ((chunk[offset] & mask[0]) == mask[0]) i++; // Ist Char ein UTF-8 zeichen?
    if ((chunk[offset] & mask[1]) == mask[1]) i++;
    if ((chunk[offset] & mask[2]) == mask[2]) i++;
    for(int j = 0; j < i; bytes[j] = chunk[offset + j++]); // Zeichen einlesen
    bytes[i] = '\0'; // Letztes Byte null setzten
    return i;
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

void writeOrderedTweets(char* TWEETS) {
    char* tweet = TWEETS;
    char buffer[linesToRead * TSIZE * 3];
    for (int i = 0; i < linesToRead; i++, tweet+=TSIZE) {
        printTweet(tweet, buffer + i * TSIZE);
    }
    
}
 

void swap(char **str1_ptr, char **str2_ptr)
{
  char *temp = *str1_ptr;
  *str1_ptr = *str2_ptr;
  *str2_ptr = temp;
}  

void compare(char* TWEETS, int i, int k, int dir) {
char* tweet1 = TWEETS + i * TSIZE;
char* tweet2 = TWEETS + k * TSIZE;
        for (int j = 6; j < TSIZE; j++) {
        if (dir == (tweet1[j] > tweet2[j])) {
            
            swap(&tweet1, &tweet2);
            
            printf("Swap\n");
            
            break;
        }
    }
}

void bitonicMerge(char* TWEETS, int start, int lines, int dir) {
    if (lines > 1) {
        int k = lines / 2;
        for (int i = start; i < start + k ; i++) 
          compare(TWEETS, i, k, dir);
        
        bitonicMerge(TWEETS, start,     k, dir);
        bitonicMerge(TWEETS, start + k, k, dir);
    }
}

void bitonicSort(char* TWEETS, int start, int lines, int dir) {
    if (lines > 1) {
        int k = lines / 2;
        bitonicSort(TWEETS, start,     k, 1 /* ASCENDING */);
        bitonicSort(TWEETS, start + k, k, 0 /* DESCENDING */);
        bitonicMerge(TWEETS, start, lines, dir);
    }
}

/*
void compare2(char* tweet1, char* tweet2, int dir) {
    for (int i = 6; i < TSIZE; i++) {
        if (dir == (tweet1[i] > tweet2[i])) {
            
        }
    }
}

int compare2(const char* ptr1, const char* ptr2) {
	int i;
	char* t1 = ptr1;
	char* t2 = ptr2;
	for (i=6; i<TSIZE; i++) {
		if (t1[i] > t2[i]) return -1;
		if (t2[i] > t1[i]) return 1;
	}
	return 0;
}


void bitonicSort2(char TWEETS[]){
  for (int k = 2; k <= linesToRead; k *= 2) //Parallel bitonic sort 
  { 
    for (int j = k / 2; j>0; j /= 2) //Bitonic merge 
    { 
      int ixj = tid ^ j; //XOR 
      if (ixj > tid) 
      { 
        if ((tid & k) == 0) // ascending descending 
        { 
          if (compare2(TWEETS[tid], TWEETS[void swap(int surface[][20], int x1, int y1, int x2, int y2) {
  int temp = surface[x1][y1];
    surface[x1][y1] = surface[x2][y2];
    surface[x2][y2] = temp;
}ixj]) == 1) 
            swap(TWEETS[tid], TWEETS[ixj]);
           
        } 
        else 
        { 
          if (compare2(TWEETS[tid], TWEETS[ixj]) == -1) 
            swap(TWEETS[tid], TWEETS[ixj]); 
        } 
      }
    } 
  } 
} */


void handle_error(const char* msg) {
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
    fprintf(stderr, "invalid month: %s\n", line);
    exit(3);
}

int countHits(const char* line, const char* key) {
    int n = strlen(key);
    int k = strlen(line) - n;
    int i;
    int hits = 0;
    for (i=0; i<k; i++, line++)
        if (*line == *key)
            if (strncmp(line, key, n) == 0)
                hits++;
    return hits;
}


void writeTweet(char* TWEETS, char* tweet, const int fn, const int ln, const int hits,
                const int month, const int day, char* line) {
    short* ptr1 = (short*) tweet;
    *ptr1 = (short) fn;
    int* ptr2 = (int*) (tweet + 2);
    *ptr2 = ln;
    *(tweet+6) = (char) hits;
    *(tweet+7) = (char) month;
    *(tweet+8) = (char) day;
    int i;
    int n = TSIZE - 9;
    for (i=strlen(line); i<n; i++) line[i] = ' '; // padding
      memcpy(tweet + 9, line, n);
    
}

char isLastProc(const int rank, const int size) {
    return rank != 0 && rank == size -1;
}

void exec(const int rank, const int size, const char* key) {

    FILE* file = fopen(FIN, "r");
    if (file == NULL) handle_error("open");
    
    
    int BUFFER_SIZE = MAX_LINE_SIZE;
    linesToRead = TNUM / size;
    int globalstart = rank * linesToRead;
    int globalend   = globalstart + linesToRead;
    
    // Der letzter Prozessor
    if (isLastProc(rank, size)) {
        globalend = TNUM;
        linesToRead = TNUM - linesToRead * (size -1) ;
    };

    
    // char Tweets[lineToRead][TSIZE]
    char* TWEETS = (char*) calloc(linesToRead, TSIZE * sizeof(char));
    if(TWEETS == NULL){
       fprintf(stderr, "Not enough memory. \n");
       exit(255);
    }
    
    char buf[BUFFER_SIZE];
    
    int lines = 0;
    char* tweet;
    char* line;
    for(tweet = TWEETS; (line = fgets(buf, BUFFER_SIZE, file)) != NULL; ++lines) {
        if(lines < globalstart) continue;
        if(lines > globalend -1) break;
        
        int fn = readNumber(&line);
        int ln = readNumber(&line);
        int month = readMonth(&line);
        int day = readNumber(&line);
        
        int hits = countHits(line, key);
        writeTweet(TWEETS, tweet, fn, ln, hits, month, day, line);
        
        tweet += TSIZE;
    }
    
    fclose(file);

    bitonicSort(TWEETS, 0, TNUM, 1);
    //bitonicSort2(TWEETS);
    writeOrderedTweets(TWEETS);
    
    free(TWEETS);
    
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
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return EXIT_SUCCESS;
}


