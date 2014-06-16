/*
 ============================================================================
 Name        : Tweet.h
 Author      : Daniil Tomilow
 ============================================================================
 */

#include "Tweet.h"



void addToken(Tweet* self, const Token* token) {
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


/**
 * Tokenize text
 */
void tokenizer(Dictionary* dict, char* tweetText, unsigned lenght, Tweet* tweet) {
    for(char *tokenStr = strtok(tweetText, " ");
            tokenStr != NULL;
            tweet->addToken(tweet, dict->put(dict, tokenStr)),
            tokenStr = strtok(NULL, " "));
}

/**
 * Tweet Beispiel:
 * {"text":"untung mbah pulang jadinya bisatanya","created_at":"Thu Feb 06 09:25:44 +0000 2014"}
 *
 * Diese Funktion parst von Beispiel Datei die so aussieht:
 * "0 2 Mar 08 @w0nderlxss No Hi's Mine"
 */
void parseTweet(Dictionary* dict, char* tweetData, unsigned int line, short int length, const int rank) {
    Tweet* tweet = newTweet();
    
    tweet->line = line;
    tweet->length = length;
    tweet->rank = rank;

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
            tokenizer(dict, tweetData + i , length, tweet);
            return;
            break;
        }
    }
}

Tweet* newTweet(){
    // Replace this with newTweet()
    Tweet *tweet = malloc(sizeof(Tweet));
    tweet->size = 0;
    tweet->addToken = addToken;
    
    return tweet;
}
