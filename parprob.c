/*
 * by Rik de Hoop & Nikolaos Trigonis
 * A Genetic Algorithm for the Partition Problem
 * AKA: parprob.c
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// #define DEBUG
#define MAX_COUNT 10000
#define numBlocks 20
#define chromlength numBlocks
#define popSize 10

// test for 2+ towers
const int towers = 2;
// gets set at initialize() 
int totalHeight;

int blocks[numBlocks];
bool generation[popSize][chromlength];

typedef struct {
  bool printChrom;
  bool printFullSum;
  bool printHeightSum;
  bool printGen;
  bool printGenIts;
} Options;

void readInput(int, char **);
int introduce(bool *);
int mutate(bool *);
int crossOver(bool *, bool *);
int heightOfTower(bool *, bool);
int heightDif(bool *);
void initialize();
void fullPrint(Options, int *, int, int, int);
void printGen();
void printChrom(int *);
void printHeightSum(int *, int);
void printFullSum(int *, int);
int rankGen(int *);
int halfCrossOver(bool *, bool *);
int towerRand();
int fitness(int *);

int main(int argc, char **argv) {
  // default options
  Options options = {1, 1, 0, 0, 1};
  int *ranks = malloc(popSize * sizeof(int));
  srand(time(NULL));
  readInput(argc, argv);
  initialize();
  
  int oldDif = -1;
  int newDif;
  int cnt = 0;
  int gen = 0;

  // main loop
  do {
    // rank the individuals and get minDif
    newDif = rankGen(ranks);
    if (oldDif == newDif) {
      cnt++;
    } else {
      cnt = 0;
      oldDif = newDif;
    }
    for (int i = popSize / 5; i < popSize; i++) {
      mutate(generation[ranks[i]]);
    }
    for (int i = popSize / 2; i < popSize; i++) {
      halfCrossOver(generation[ranks[i]], generation[ranks[0]]);
    }
#ifdef DEBUG
    printf("gen %d: %d\n", gen, newDif);
#endif
    gen++;

  } while ((cnt < MAX_COUNT) && (newDif != 0));

  fullPrint(options, ranks, newDif, gen, cnt);
  free(ranks);
  return 0;
}
// returns the fitness of a chromosome, for every doubling of towers, the max deviation increases by 50%
int fitness(int *chromosome) {
  int optimum = totalHeight/towers;
  int deviation = 0, height = 0;
  // find the deviation of the optimum for every tower
  for (int t = 0; t < towers; t++) {
    for (int g = 0; g < chromlength; g++) {
      if (chromosome[g] == t) {
        height += blocks[g];
      }
    }
    deviation += abs(height-optimum);
    height = 0;
  }
  return deviation;
}
// picks a random tower, divide RAND_MAX by towers and divide
// rand() by that to get a random number, as this method has a low chance of
// getting towers, in case this happens the tower index decrements
int towerRand() {
  int tower = rand() / (RAND_MAX / towers);
  if (tower == towers) {
      tower--;
  }
  return tower;
}
// print final results
void fullPrint(Options options, int *ranks, int minDif, int gen, int cnt) {
  if (options.printChrom) {
    printChrom(ranks);
  }
  if (options.printFullSum) {
    printFullSum(ranks, minDif);
  }
  if (options.printHeightSum) {
    printHeightSum(ranks, minDif);
  }
  if (options.printGen) {
    printGen();
  }
  if (options.printGenIts) {
    printf("Generation: %d; Iteration: %d\n", gen - cnt - 1, gen - 1);
  }
}

// returns the lowest difference between heights and ranks the individuals
int rankGen(int *ranks) {
  int temp;
  int *rankedDifs = malloc(popSize * sizeof(int));
  for (int i = 0; i < popSize; i++) {
    ranks[i] = i;
  }
  for (int i = 0; i < popSize; i++) {
    rankedDifs[i] = heightDif(generation[i]);
  }
  for (int i = 0; i < popSize - 1; i++) {
    for (int j = i + 1; j < popSize; j++) {
      if (rankedDifs[i] > rankedDifs[j]) {
        temp = rankedDifs[i];
        rankedDifs[i] = rankedDifs[j];
        rankedDifs[j] = temp;

        temp = ranks[i];
        ranks[i] = ranks[j];
        ranks[j] = temp;
      }
    }
  }
  temp = rankedDifs[0];
  free(rankedDifs);
  return temp;
}

void printGen() {
  for (int i = 0; i < popSize; i++) {
    for (int j = 0; j < chromlength; j++) {
      printf("%d ", generation[i][j]);
    }
    printf("\n");
  }
}

// fills all chromosomes of the first generation
void initialize() {
  int total = 0;
  for (int i = 0; i < numBlocks; i++) {
    total += blocks[i];
  }
  totalHeight = total;
  for (int i = 0; i < popSize; i++) {
    introduce(generation[i]);
  }
}

// reads heights or generates them
void readInput(int argc, char **argv) {
  // implement something that reads a file later, and a randomly filled one
  if (argc == 1) {
    for (int i = 0; i < numBlocks; i++) {
      (void)scanf("%d ", &blocks[i]);
    }
  } else if (!strcmp(argv[1], "rand")) {
    // random integers
    printf("testestrandrandrand\n");
  }
}

// fill a chromosome with random values
int introduce(bool *chromosome) {
  for (int g = 0; g < chromlength; g++) {
    chromosome[g] = towerRand();
  }
}

// mutates 1 randome gene in a chromosome
int mutate(bool *chromosome) {
  int idx = rand() % chromlength;
  chromosome[idx] = towerRand();
}

// after a random index, swap the values of both chromosomes
int crossOver(bool *chrom1, bool *chrom2) {
  bool temp;
  int start = (rand() % (chromlength - 1)) + 1;
  for (int i = start; i < chromlength; i++) {
    temp = chrom1[i];
    chrom1[i] = chrom2[i];
    chrom2[i] = temp;
  }
}

int halfCrossOver(bool *target, bool *origin) {
  int start = (rand() % (chromlength - 1)) + 1;
  for (int i = start; i < chromlength; i++) {
    target[i] = origin[i];
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

void printFullSum(int *ranks, int minDif) {
  int first = 1;
  printf("T set:\n");
  for (int i = 0; i < chromlength; i++) {
    if (generation[ranks[0]][i]) {
      if (first) {
        first = 0;
        printf("%d", blocks[i]);
      } else {
        printf("+%d", blocks[i]);
      }
    }
  }
  printf("=%d\nF set:\n", heightOfTower(generation[ranks[0]], true));
  for (int i = 0; i < chromlength; i++) {
    if (!generation[ranks[0]][i]) {
      printf("-%d", blocks[i]);
    }
  }
  printf("=-%d\nabs(%d - %d) = %d\n",
         heightOfTower(generation[ranks[0]], false),
         heightOfTower(generation[ranks[0]], true),
         heightOfTower(generation[ranks[0]], false), minDif);
}

void printChrom(int *ranks) {
  printf("Chromosome: ");
  for (int i = 0; i < chromlength; i++) {
    printf("%d", generation[ranks[0]][i]);
  }
  printf("\n");
}

void printHeightSum(int *ranks, int minDif) {
  int trueTower = heightOfTower(generation[ranks[0]], true);
  int falseTower = heightOfTower(generation[ranks[0]], false);
  printf("%d - %d = ", (trueTower > falseTower) ? trueTower : falseTower,
         (trueTower > falseTower) ? falseTower : trueTower);
  printf("%d\n", minDif);
}