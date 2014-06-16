/*
 ============================================================================
 Name        : Dictionary.h
 Author      : Daniil Tomilow
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ENTRIES 500 // Max eintraege im Dictionary

typedef struct Token {
    char* value; // der eigentliche token (komprimieren?)
    unsigned int length; // Token laenge
    unsigned int used; // Anzahl wie oft das Token benutzt wird
} Token;


typedef struct Dictionary {
    unsigned int size; // Anzahl der Woerter im Dict
    const Token* (*put)(struct Dictionary*, char*);
    void (*printDict)(const struct Dictionary*);
    Token* tokens[ENTRIES];
} Dictionary;

/**
 * Constructor: Returns new Dictionary
 */
Dictionary* newDictionary();

//void destroyDictionary(Dictionary *dict);

/**
 * Add new unique Token to Dictionary.
 * Returns pointer to new or existing Token.
 */ 
const Token* put(Dictionary *self, char *value);

/**
 * Prints the Dictionary.
 */
void printDict(const Dictionary * dict);

