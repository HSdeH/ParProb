/*
 * by Rik de Hoop & Nikolaos Trigonis
 * A Genetic Algorithm for the Partition Problem
 * AKA: parprob.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>

#define numBlocks 20
#define chromlength numBlocks
#define popSize 10

int blocks[numBlocks];
bool generation[popSize][chromlength];

void readHeights(int, char **);
int fillChromosome(bool *);
int mutate(bool *);
int crossOver(bool *, bool *);
int heightOfTower(bool *, bool);
int heightDif(bool *);

int main(int argc, char **argv) {
  srand(time(NULL));  
  readHeights(argc, argv);
  
  // fitness, necessary or not, just inverse of height dif
  // are you allowed to initialize with (a) heuristic(s)? 
  return 0;
}

// reads heights or generates them
void readHeights(int argc, char **argv) {
  // implement something that reads a file later, and a randomly filled one
  if (argc == 1) {
    printf("Insert %d numbers:\n", numBlocks);
    for (int i = 0; i < numBlocks; i++) {
      (void)scanf("%d ", &blocks[i]);
    }
  } else if (!strcmp(argv[1], "rand")) {
    // random integers
    printf("testestrandrandrand\n");
  } else {

  }
}

// fill a chromosome with random values
int fillChromosome(bool *chromosome) {
  for (int g = 0; g < chromlength; g++) {
    chromosome[g] = (rand() > RAND_MAX / 2);
  }
}

// flips 1 randome gene in a chromosome
int mutate(bool *chromosome) {
  int idx = rand() % chromlength;
  chromosome[idx] = !chromosome[idx]; 
}

// after a random index, swap the values of both chromosomes
int crossOver(bool *chrom1, bool *chrom2) {
  bool temp;
  int start = (rand() % (chromlength-1)) + 1;
  for (int i = start; i < chromlength; i++) {
    temp = chrom1[i];
    chrom1[i] = chrom2[i];
    chrom2[i] = temp;
  }
}

// returns the height of T or F within a chromosome
int heightOfTower(bool *chromosome, bool b) {
  int height = 0;
  for (int i = 0; i < chromlength; i++) {
    if (chromosome[i] == b) {
      height += blocks[i];
    }
  }
  return height;
}

// returns the height difference of a chromosome
int heightDif(bool *chromosome) {
  int dif = 0;
  for (int i = 0; i < chromlength; i++) {
    if (chromosome[i]) {
      dif += blocks[i];
    } else {
      dif -= blocks[i];
    }
  }
  return abs(dif);
}
