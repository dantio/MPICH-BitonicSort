/*
 ============================================================================
 Name        : Dictionary.c
 Author      : Daniil Tomilow
 ============================================================================
 */

#include "Dictionary.h"

/**
 * Add new unique Token to Dictionary.
 * Returns pointer to new or existing Token.
 */
const Token* put(Dictionary *self, char *value) {
    if(self->size + 1 > ENTRIES) {
        printf("Use realloc here\n");
        return NULL;
    }
    // 1. suche nach den gleichen String
    unsigned int length = strlen(value);
    if(self->size > 0) {
        for(int i = 0; i < self->size; i++) {
            Token* comT = self->tokens[i];
            if(comT->length != length && strncmp(comT->value, value, length) != 0)
                continue;
            else {
                // Gleicher Token gefunden!
                comT->used++;
                return comT;
            }
        }
    }
    // 2. Kein gleicher Token, erstelle einen neuen
    Token* token = newToken(value);
    
    self->tokens[self->size] = token;
    self->size++;
    return token;
}

void printDict(const Dictionary * dict) {
    if(dict == NULL) return;
    printf("== Rank: %d; Size: %d; Tokens:\n", dict->rank,  dict->size);
    if(dict->size > 0) {
        Token * t;
        for(int i = 0; i < dict->size; i++) {
            t = dict->tokens[i];
           // if(t->used > 1)
                printf("Used: %d; Value: \"%s\" \n", t->used, t->value);
        }
    }
}

Dictionary* newDictionary(int rank) {
    Dictionary* dict = malloc(sizeof(Dictionary));
    if(dict == NULL) {
          printf("No more memory \n");
          return NULL;
    }
    dict->rank = rank;
    // Assign methods
    dict->put = put;
    dict->printDict = printDict;
    
    dict->size = 0;
    
    return dict;
};

void destroyDictionary(Dictionary *dict) {
    // Todo
    
}

Token* newToken(char* value) {
    unsigned int length = strlen(value);
    Token* newToken = malloc(sizeof(Token));
    newToken->value = malloc(sizeof(char) * length);
    strncpy(newToken->value, value, length);
    newToken->length = length;
    newToken->used = 1;
    
    return newToken;
}


