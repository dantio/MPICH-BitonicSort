/*
 ============================================================================
 Name        : Tweetonic.c
 Author      : Daniil Tomilow
 ============================================================================
 */
#include "Dictionary.h"
#include "Tweet.h"

#include <mpi.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_SIZE 652 // 652 / 4 = 163

// mask values for bit pattern of first byte in multi-byte
// UTF-8 sequences:
// 192 - 110xxxxx - for U+0080 to U+07FF
// 224 - 1110xxxx - for U+0800 to U+FFFF
// 240 - 11110xxx - for U+010000 to U+1FFFFF
#define U_MAX_BYTES 4
static unsigned short mask[] = {192, 224, 240};

unsigned int lineOffset = 0;
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


void exec(MPI_File* in, const int rank, const int size, const int overlap) {
    Dictionary* dict = newDictionary();
    
    MPI_Offset globalstart;
    int lastProc = size - 1; // Letzter Prozess
    int mysize;
    char* chunk;
    {
        MPI_Offset globalend;
        MPI_Offset filesize;
        // Datei groese ermitteln
        MPI_File_get_size(*in, &filesize);
        // Letztes null byte ignorieren
        filesize--;
        // Jeder Prozessor bekommt ein Teil der Datei
        mysize = filesize / size;
        globalstart = rank * mysize;
        globalend   = globalstart + mysize - 1;
        // Der letzter Prozessor
        if (rank == lastProc) globalend = filesize - 1;
        // Ueberlappen um keine halben Zeilen zu haben
        // Beim letzten nicht, weil er bis zum Ende der Datei liest
        if (rank != lastProc) globalend += overlap;
        mysize =  globalend - globalstart + 1;
        // Speicher alloziieren
        chunk = malloc((mysize + 1) * sizeof(char));
        // Chunks einlesen
        MPI_File_read_at_all(
            *in,
            globalstart,
            chunk,
            mysize,
            MPI_CHAR,
            MPI_STATUS_IGNORE);
        chunk[mysize] = '\0';
    }
    //
    // Suche den "richtigen" Start. (Weil es moeglich ist mitten in der Zeile
    // anzufangen oder aufzuhoeren).
    //
    int locstart = 0; // Von
    int locend = mysize - 1; // Bis
    // Rank 0 liest von anfang der Datei, der Start muss nicht gesucht werde
    if (rank != 0) {
        while(chunk[locstart] != '\n') locstart++;
        locstart++;
    }
    // Der letzte liest einfach bis zum ende
    if (rank != lastProc) {
        locend -= overlap;
        while(chunk[locend] != '\n') locend++;
    }
    mysize = locend - locstart + 1;
    // Zeile lesen, parsen
    // pthread_t thread1;
    char tweet[MAX_LINE_SIZE];
    unsigned int j = 0;
    unsigned int lineNum = 0;
    char c[U_MAX_BYTES];
    
    // Hier koennte openMP stehen
    for (unsigned int i = locstart; i <= locend;) {
        char* cp = &c;
        unsigned int charLen = u_getc(chunk, i, c);
        // Addiere die Laenge des UTF-8 chars
        i += charLen;
        // Kopiere UTF-8 Zeichen so lange bis Zeichen zu ende ist
        if(c[charLen - 1] != '\n' && i + 1 < locend) while(*cp) tweet[j++] = *cp++;
        // Zeile ist zu ende
        else {
            tweet[j] = '\0';
            parseTweet(dict, tweet, lineNum, j, rank);
            j = 0;
            lineNum++;
            //if(thread1 != NULL) pthread_join(thread1, NULL);
            //pthread_create(&thread1, NULL, (void *) &parseTweet, (char *) &tweet);
        }
    }
    // Zeilennummer austauschen
    if(rank != lastProc) MPI_Send(&lineNum, 1, MPI_UNSIGNED, rank + 1, 0 /* TAG */, MPI_COMM_WORLD);
    if(rank != 0) MPI_Recv(&lineOffset, 1, MPI_UNSIGNED, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    dict->printDict(dict);
}



int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(NULL, NULL);
    
    // Get the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_File fh;
    if(MPI_SUCCESS != MPI_File_open(
                MPI_COMM_WORLD,
                "twitter.data10",
                MPI_MODE_RDONLY,
                MPI_INFO_NULL, &fh) ) {
        fprintf(stderr, "Couldn't open file ");
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    exec(&fh, rank, size, MAX_LINE_SIZE);
    MPI_File_close(&fh);
    // Print off a hello world message
//  printf("Hello world from processor %s, rank %d out of %d processors\n",
    //       processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return EXIT_SUCCESS;
}


