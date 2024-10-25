/*
 * by Rik de Hoop & Nikolaos Trigonis
 * A Genetic Algorithm for the Partition Problem
 * AKA: parprob.c
 */

#include <math.h>
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

int towerAmount = 2;
int totalHeight;

int blocks[numBlocks];
int generation[popSize][chromlength];

typedef struct {
  bool printChrom;
  bool printFullSum;
  bool printHeightSum;
  bool printGen;
  bool printGenIts;
} Options;

char toBase64(int);
void readInput(int, char **, Options *);
void introduce(int *);
void mutate(int *);
void crossOver(int *, int *);
int heightOfTower(int *, int);
int heightDif(int *);
void initialize();
void fullPrint(Options *, int *, double, int, int);
void printGen();
void printChrom(int *);
void printHeightSum(int *, int);
void printFullSum(int *, double);
double rankGen(int *);
void halfCrossOver(int *, int *);
int towerRand();
double deviation(int *);

int main(int argc, char **argv) {
  Options options = {1, 1, 0, 0, 1}; // default options
  int *ranks = malloc(popSize * sizeof(int));
  srand(time(NULL));
  readInput(argc, argv, &options);
  initialize();

  double oldDev = -1;
  double newDev;
  int cnt = 0;
  int gen = 0;

  // main loop
  do {
    // rank the individuals and get the lowest deviation
    newDev = rankGen(ranks);
    if (oldDev == newDev) {
      cnt++;
    } else {
      cnt = 0;
      oldDev = newDev;
    }
    // change the new generation
    for (int i = popSize / 5; i < popSize; i++) {
      mutate(generation[ranks[i]]);
    }
    for (int i = popSize / 2; i < popSize; i++) {
      halfCrossOver(generation[ranks[i]], generation[ranks[0]]);
    }
    gen++;
  } while ((cnt < MAX_COUNT) &&
           (newDev != 0));  // continue until a solution has been found or
                            // nothing has been found for a while

  fullPrint(&options, ranks, newDev, gen, cnt);
  free(ranks);
  return 0;
}
// implemented to make the printed chromosome better capable of handling lots of
// towers bit too simple right now
char toBase64(int n) {
  char base64[] = {
      "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+/"};
  return base64[n % 64];
}
// returns the average deviation of a chromosome, for every doubling of towers,
// the lower the deviation the higher the fitness
double deviation(int *chromosome) {
  int *towers = calloc(sizeof(int), towerAmount);
  double optimum = (double)totalHeight / (double)towerAmount;
  double deviation = 0.0;
  // build every tower
  for (int g = 0; g < chromlength; g++) {
    towers[chromosome[g]] += blocks[g];
  }
  // then add up the deviation from the optimum for every tower
  for (int t = 0; t < towerAmount; t++) {
    deviation += fabs((double)towers[t] - optimum);
  }
  free(towers);
  return deviation / towerAmount;
}
// picks a random tower, divide RAND_MAX by towers and divide
// rand() by that to get a random number, as this method has a low chance of
// getting towers, in case this happens the tower index decrements
int towerRand() {
  int tower = rand() / (RAND_MAX / towerAmount);
  return tower == towerAmount ? tower - 1 : tower;
}
// print final results
void fullPrint(Options *options, int *ranks, double deviation, int gen,
               int cnt) {
  if (options->printChrom) {
    printChrom(ranks);
  }
  if (options->printFullSum) {
    printFullSum(ranks, deviation);
  }
  if (options->printHeightSum) {
    printHeightSum(ranks, deviation);
  }
  if (options->printGen) {
    printGen();
  }
  if (options->printGenIts) {
    printf("Generation: %d; Iteration: %d\n", gen - cnt - 1, gen - 1);
  }
}

// returns the lowest difference between heights and ranks the individuals
double rankGen(int *ranks) {
  double temp;
  double *rankedDifs = malloc(popSize * sizeof(double));
  for (int i = 0; i < popSize; i++) {
    ranks[i] = i;
  }
  for (int i = 0; i < popSize; i++) {
    rankedDifs[i] = deviation(generation[i]);
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
      printf("%c", toBase64(generation[i][j]));
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
void readInput(int argc, char **argv, Options *options) {
  // implement something that reads a file later, and a randomly filled one
  int randarg = 0;
  int a = 1;
  while (a != argc) {
    if (!strcmp(argv[a], "-pg")) {
      options->printGen = true;
    }
    if (!strcmp(argv[a], "rand")) {
      // instead of determining the output yourself, let every block have a
      // height of min-max
      randarg = 1;
      int randMin = atoi(argv[a + 1]), randMax = atoi(argv[a + 2]) + 1;

      for (int i = 0; i < numBlocks; i++) {
        blocks[i] = (rand() % (randMax - randMin)) + randMin;
      }
      a = a + 2;
    }
    a++;
  }
  if (!randarg) {
    for (int i = 0; i < numBlocks; i++) {
      (void)scanf("%d ", &blocks[i]);
    }
  }
}

// fill a chromosome with random values
void introduce(int *chromosome) {
  for (int g = 0; g < chromlength; g++) {
    chromosome[g] = towerRand();
  }
}

// mutates 1 randome gene in a chromosome
void mutate(int *chromosome) {
  int idx = rand() % chromlength;
  int rndm;
  do {
    rndm = towerRand();
  } while (rndm == chromosome[idx]);
  chromosome[idx] = rndm;
}

// after a random index, swap the values of both chromosomes
void crossOver(int *chrom1, int *chrom2) {
  bool temp;
  int start = (rand() % (chromlength - 1)) + 1;
  for (int i = start; i < chromlength; i++) {
    temp = chrom1[i];
    chrom1[i] = chrom2[i];
    chrom2[i] = temp;
  }
}
// crosOver() but the origin doesn't get changed
void halfCrossOver(int *target, int *origin) {
  int start = (rand() % (chromlength - 1)) + 1;
  for (int i = start; i < chromlength; i++) {
    target[i] = origin[i];
  }
}

// returns the height of T or F within a chromosome
int heightOfTower(int *chromosome, int t) {
  int height = 0;
  for (int i = 0; i < chromlength; i++) {
    if (chromosome[i] == t) {
      height += blocks[i];
    }
  }
  return height;
}

// returns the height difference of a chromosome, assumes 2 towers
int heightDif(int *chromosome) {
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

void printFullSum(int *ranks, double deviation) {
  int first, sum;
  for (int s = 0; s < towerAmount; s++) {
    sum = 0;
    first = 1;
    printf("Set %c:\t", toBase64(s));
    for (int b = 0; b < chromlength; b++) {
      if (generation[ranks[0]][b] == s) {
        if (first) {
          first = 0;
          printf("%d", blocks[b]);
        } else {
          printf("+%d", blocks[b]);
        }
        sum += blocks[b];
      }
    }
    printf("=%d; dev: %.3f\n", sum,
           (double)sum - ((double)totalHeight / (double)towerAmount));
  }
  printf("Optimal: %.3f; Average deviation: %.3f; Total deviation: %.3f\n",
         (double)totalHeight / (double)towerAmount, deviation,
         deviation * towerAmount);
}

void printChrom(int *ranks) {
  printf("Chromosome: ");
  for (int i = 0; i < chromlength; i++) {
    printf("%c", toBase64(generation[ranks[0]][i]));
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