/*
 ============================================================================
 Name        : Tweetonic.c
 Author      : Eric Wündisch, Daniil Tomilow
 ============================================================================
 */
#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_SIZE 652 // 652 / 4 = 163
#define TOKENS 75 // Max tokens in einem Tweet
#define ENTRIES 500 // Max eintraege im Dictionary

// see http://www.mediabistro.com/alltwitter/files/2010/04/twitter_know_your_limits.png
#define MAX_TWEET_TEXT 140
#define MAX_USER_NAME 15

// mask values for bit pattern of first byte in multi-byte
// UTF-8 sequences:
// 192 - 110xxxxx - for U+0080 to U+07FF
// 224 - 1110xxxx - for U+0800 to U+FFFF
// 240 - 11110xxx - for U+010000 to U+1FFFFF
#define U_MAX_BYTES 4
static unsigned short mask[] = {192, 224, 240};

typedef struct Token {
    char value[MAX_TWEET_TEXT * U_MAX_BYTES]; // der eigentliche token (komprimieren?)
    unsigned int length; // Token laenge
    unsigned int used; // Anzahl wie oft das Token benutzt wird
} Token;


typedef struct Dictionary {
    unsigned int size; // Anzahl der Woerter im Dict
    const Token * (*put)(struct Dictionary *, char *);
    Token *tokens[ENTRIES];
} Dictionary;

typedef struct Tweet {
    unsigned int line; // Zeilennummer/ ID in der Twitter Datei
    unsigned int date; // Datum, tagesgenau
    unsigned int length; // Tweet laenge
    unsigned int size; // Anzahl der Tokens
    int rank; // Rang on Prozessoer der diesen Tweet ausgewertet hat
    char createdBy[MAX_USER_NAME]; // Erstellt von
    
    void (*addToken)(struct Tweet *, const Token *);
    Token *tokens[TOKENS]; // Zeiger auf tokens id
} Tweet;


// Alle tweets im Chunk von Prozessor
Tweet *tweets[500];

// Dict funktionen
// Woerterbuch
Dictionary dict;
const Token * put(Dictionary *self, char *value) {
    if(self->size + 1 > ENTRIES) {
        printf("Use realloc here\n");
        return;
    }
    // 1. suche nach den gleichen String
    unsigned int length = strlen(value);
    if(self->size > 0) {
        for(int i = 0; i < self->size; i++) {
            Token * comT = self->tokens[i];
            if(comT->length != length && strncmp(comT->value, value, length) != 0)
                continue;
            else {// Gleicher Token gefunden!
                comT->used++;
                return comT;
            }
        }
    }
    // 2. Kein gleicher Token, erstelle einen neuen
    Token *newToken = malloc(sizeof(Token));
    strcpy(&newToken->value, value);
    newToken->length = length;
    newToken->used = 1;
    self->tokens[self->size] = newToken;
    self->size++;
    return newToken;
}

void printDict(Dictionary * dict) {
    if(dict == NULL) return;
    printf("Size: %d; Tokens:\n", dict->size);
    if(dict->size > 0) {
        Token * t;
        for(int i = 0; i < dict->size; i++) {
            t = dict->tokens[i];
            if(t->used > 100)
            printf("Used: %d; Value: \"%s\" \n", t->used, t->value);
        }
    }
}

// Tweet funktionen
void addToken(Tweet *self, const Token * token) {
    if(self->size + 1 > TOKENS) {
        printf("Too much tokens.\n");
        return;
    }
    self->tokens[self->size] = token;
    self->size++;
}


void printTweet(Tweet * tweet) {
    if(tweet == NULL) return;
    printf("Line: %d; Rank: %d; Tweet: {\"", tweet->line, tweet->rank);
    if(tweet->size > 0) {
        Token * t;
        for(int i = 0; i < tweet->size; i++) {
            t = tweet->tokens[i];
            printf("%s ", t->value);
        }
        printf("\"} \n");
    }
}



unsigned int lineOffset = 0;
/**
 * UTF-8. Return length of the UTF-8 Character.
 */
int u_getc(char *chunk, int offset, char *bytes) {
    int i = 1;
    if ((chunk[offset] & mask[0]) == mask[0]) i++; // Ist Char ein UTF-8 zeichen?
    if ((chunk[offset] & mask[1]) == mask[1]) i++;
    if ((chunk[offset] & mask[2]) == mask[2]) i++;
    for(int j = 0; j < i; bytes[j] = chunk[offset + j++]); // Zeichen einlesen
    bytes[i] = '\0'; // Letztes Byte null setzten
    return i;
}


/**
 * Tokenize text
 */
void tokenizer(char *tweetText, unsigned lenght, Tweet *tweet) {
    for(char *tokenStr = strtok(tweetText, " ");
            tokenStr != NULL;
            tweet->addToken(tweet, dict.put(&dict, tokenStr)),
            tokenStr = strtok(NULL, " "));
}

/**
 * Tweet Beispiel:
 * {"text":"untung mbah pulang jadinya bisatanya","created_at":"Thu Feb 06 09:25:44 +0000 2014"}
 *
 * Diese Funktion parst von Beispiel Datei die so aussieht:
 * "0 2 Mar 08 @w0nderlxss No Hi's Mine"
 */
void parseTweet(char *tweetData, unsigned int line, short int length, const int rank) {
    // Replace this with newTweet()
    Tweet *tweet = malloc(sizeof(Tweet));
    tweet->line = line;
    tweet->length = length;
    tweet->rank = rank;
    tweet->size = 0;
    tweet->addToken = addToken;
    int pos = 0; // position von whitespace
    for(int i = 0; i < length; i++) {
        if(tweetData[i] == ' ') {
            pos++;
            continue;
        }
        // Tweettext
        switch(pos) {
        case 0:
            break; // Ignore
        case 1:
            break; // ID/Line number
        case 2:
            break; // Monat
        case 3:
            break; // Tag
        case 4: // Tokenizer
            tokenizer(tweetData + i , length, tweet);
            //printTweet(tweet);
            return;
            break;
        }
    }
}

void exec(MPI_File *in, const int rank, const int size, const int overlap) {
    MPI_Offset globalstart;
    int lastProc = size - 1; // Letzter Prozess
    int mysize;
    char *chunk;
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
        char *cp = &c;
        unsigned int charLen = u_getc(chunk, i, &c);
        // Addiere die Laenge des UTF-8 chars
        i += charLen;
        // Kopiere UTF-8 Zeichen so lange bis Zeichen zu ende ist
        if(c[charLen - 1] != '\n' && i + 1 < locend) while(*cp) tweet[j++] = *cp++;
        // Zeile ist zu ende
        else {
            tweet[j] = '\0';
            parseTweet(tweet, lineNum, j, rank);
            j = 0;
            lineNum++;
            //if(thread1 != NULL) pthread_join(thread1, NULL);
            //pthread_create(&thread1, NULL, (void *) &parseTweet, (char *) &tweet);
        }
    }
    // Zeilennummer austauschen
    if(rank != lastProc) MPI_Send(&lineNum, 1, MPI_UNSIGNED, rank + 1, 0 /* TAG */, MPI_COMM_WORLD);
    if(rank != 0) MPI_Recv(&lineOffset, 1, MPI_UNSIGNED, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printDict(&dict);
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
    // Replace this with newDict
    dict.put = put;
    dict.size = 0;
    exec(&fh, rank, size, MAX_LINE_SIZE);
    MPI_File_close(&fh);
    // Print off a hello world message
//  printf("Hello world from processor %s, rank %d out of %d processors\n",
    //       processor_name, rank, size);
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return EXIT_SUCCESS;
}

