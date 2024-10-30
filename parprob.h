#ifndef PARPROB_C
#define PARPROB_C

#include <time.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct {
   bool defaultOptions;
   bool printChrom;
   bool printFullSum;
   bool printHeightSum;
   bool printGen;
   bool printGenIts;
   bool printTime;
   bool improvement;
   bool greedy;
} Options;

char toBase64(int);
void readInput(int, char **, Options *);
void introduce(int *);
void mutate(int *, double);
void crossOver(int *, int *);
void initialize();
void fullPrint(Options *, int *, double, int, int, int);
void printGen();
void printChrom(int *);
void printFullSum(int *, double);
double rankGen(int *);
int towerRand();
double deviation(int *);
void printTime(clock_t);
void greedy(int *);
void recombine(int *, int *);
void readFile(FILE *);

#endif