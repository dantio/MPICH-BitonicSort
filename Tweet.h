/*
 ============================================================================
 Name        : Tweet.h
 Author      : Daniil Tomilow
 ============================================================================
 */
 
#include "Dictionary.h"
#include <stdio.h>

 // see http://www.mediabistro.com/alltwitter/files/2010/04/twitter_know_your_limits.png
#define MAX_TWEET_TEXT 140
#define MAX_USER_NAME 15
#define TOKENS 75 // Max tokens in einem Tweet

typedef struct Tweet {
    unsigned int line; // Zeilennummer/ ID in der Twitter Datei
    unsigned int date; // Datum, tagesgenau
    unsigned int length; // Tweet laenge
    unsigned int size; // Anzahl der Tokens
    int rank; // Rang on Prozessoer der diesen Tweet ausgewertet hat
    char createdBy[MAX_USER_NAME]; // Erstellt von
    
    void (*addToken)(struct Tweet *, const Token*);
    Token* tokens[TOKENS]; // Zeiger auf tokens id
} Tweet;

Tweet* newTweet();

void tokenizer(Dictionary* dict, char* tweetText, unsigned lenght, Tweet* tweet);

void parseTweet(Dictionary* dict, char* tweetData, unsigned int line, short int length, const int rank);
