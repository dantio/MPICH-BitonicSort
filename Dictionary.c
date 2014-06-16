/*
 ============================================================================
 Name        : Dictionary.c
 Author      : Daniil Tomilow
 ============================================================================
 */
 
#include "Dictionary.h"

Dictionary* newDictionary(){
    Dictionary* dict = malloc(sizeof(Dictionary));
    // Assign methods
    dict->put = put;
    dict->printDict = printDict;
    dict->size = 0;
    
    return dict;
};

const Token* put(Dictionary *self, char *value) {
    if(self->size + 1 > ENTRIES) {
        printf("Use realloc here\n");
        return NULL;
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
    Token* newToken = malloc(sizeof(Token));
    newToken->value = malloc(sizeof(char) * length);
    strcpy(newToken->value, value);
    newToken->length = length;
    newToken->used = 1;
    self->tokens[self->size] = newToken;
    self->size++;
    return newToken;
}

void printDict(const Dictionary * dict) {
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
