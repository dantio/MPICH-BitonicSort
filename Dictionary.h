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

#ifndef TOKEN
#define TOKEN
typedef struct Token {
    char* value; // der eigentliche token (komprimieren?)
    unsigned int length; // Token laenge
    unsigned int used; // Anzahl wie oft das Token benutzt wird
} Token;
#endif

#ifndef DICT
#define DICT
typedef struct Dictionary {
    unsigned int size; // Anzahl der Woerter im Dict
    const Token* (*put)(struct Dictionary*, char*);
    void (*printDict)(const struct Dictionary*);
    Token* tokens[ENTRIES];
} Dictionary;
#endif

/**
 * Dictionary constructor: Returns new Dictionary
 */
Dictionary* newDictionary();

/**
 * Dictionary deconstructor: delete dictionary and every Token
 */
void destroyDictionary(Dictionary *dict);

/**
 * Token constructor: Returns new Token
 */
Token* newToken(char* value);

