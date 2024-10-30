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
int popSize = 10;
int maxDrift = 10000;
int subsets;
int numBlocks, chromlength;
int totalHeight;
int *blocks;
int **generation;

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
int heightOfTower(int *, int);
void initialize();
void fullPrint(Options *, int *, double, int, int, int);
void printGen();
void printChrom(int *);
void printFullSum(int *, double);
double rankGen(int *);
void halfCrossOver(int *, int *);
int towerRand();
double deviation(int *);
void printTime(clock_t);
void greedy(int *);
void recombine(int *, int *);

int main(int argc, char **argv) {
   Options options = {1};  // options, mostly for printing results
   srand(time(NULL));
   readInput(argc, argv, &options);
   initialize();
   if (options.greedy) {
      greedy(generation[0]);
   }
   int *ranks = malloc(popSize * sizeof(int));
   double oldDev = -1, newDev, mutationFactor;
   int gen = 0, drift = 0;
   clock_t start = clock(), end;
   // main loop
   do {
      // rank the individuals and get the lowest deviation
      newDev = rankGen(ranks);
      if (oldDev == newDev) {
         drift++;
      } else {
         if (options.improvement || options.defaultOptions) {
            printf("Gen %-7d: %.3f\n", gen,
                   newDev);  // each time an improvement has been found, print
                             // the generation and current improvement
         }
         drift = 0;
         oldDev = newDev;
      }
      mutationFactor = (double)drift / (double)maxDrift;
      // change the new generation
      for (int i = 1; i < popSize; i++) {
         recombine(generation[ranks[i]], generation[ranks[0]]);
         mutate(generation[ranks[i]], 0.2*mutationFactor); // this allows the mutation to increase to 25% of genes as drift increases
      }
      gen++;
   } while ((drift < maxDrift) &&
            (newDev != 0));  // continue until a solution has been found or
                             // maxDrift has been reached
   end = clock();
   // print everything you want to see
   fullPrint(&options, ranks, newDev, gen, drift, end - start);
   // free all memory
   free(ranks);
   free(blocks);
   for (int i = 0; i < popSize; i++) {
      free(generation[i]);
   }
   free(generation);
   return 0;
}
// implemented to make the printed chromosome better capable of handling lots of
// towers, when larger, starts from 0 again, which I do not think is a big deal
char toBase64(int n) {
   char base64[] = {
       "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+/"};
   return base64[n % 64];
}

int compare(const void *a, const void *b) { return (*(int *)b - *(int *)a); }
void greedy(int *chrom) {
   int *sums = calloc(sizeof(int), subsets);
   qsort(blocks, numBlocks, sizeof(int), compare);

   for (int g = 0; g < chromlength; g++) {
         chrom[g] = 0;
         for (int j = 1; j < subsets; j++) {
            if (sums[chrom[g]] > sums[j]) {
               chrom[g] = j;
            }
         }
         sums[chrom[g]] += blocks[g];
      
   }
   free(sums);
}

// a uniform crossOver;
void recombine(int *target, int *origin) {
  for (int g = 0; g < chromlength; g++) {
     if (rand() > RAND_MAX / 2) {
      target[g] = origin[g];
     }
  }
}

// returns the average deviation of a chromosome, for every doubling of towers,
// the lower the deviation the higher the fitness
// alternate: largest tower - optimum
double deviation(int *chromosome) {
   int *towers = calloc(sizeof(int), subsets);
   double optimum = (double)totalHeight / (double)subsets;
   double deviation = 0.0;
   // build every tower
   for (int g = 0; g < chromlength; g++) {
      towers[chromosome[g]] += blocks[g];
   }
   // then add up the deviation from the optimum for every tower
   for (int t = 0; t < subsets; t++) {
      deviation += fabs((double)towers[t] - optimum);
   }
   free(towers);
   return deviation / subsets;
}
// picks a random tower, divide RAND_MAX by towers and divide
// rand() by that to get a random number, as this method has a low chance of
// getting towers, in case this happens the tower index decrements
int towerRand() {
   int tower = rand() / (RAND_MAX / subsets);
   return tower == subsets ? tower - 1 : tower;
}
// print final results
void fullPrint(Options *options, int *ranks, double deviation, int gen, int cnt,
               int time) {
   if (options->printChrom || options->defaultOptions) {
      printChrom(ranks);
   }
   if (options->printFullSum || options->defaultOptions) {
      printFullSum(ranks, deviation);
   }
   if (options->printGen) {
      printGen();
   }
   if (options->printGenIts || options->defaultOptions) {
      printf("Generation: %d; Iteration: %d\n", gen - cnt - 1, gen - 1);
   }
   if (options->printTime || options->defaultOptions) {
      printTime(time);
   }
}

void printTime(clock_t time) {
   int s = ((double)(time) / CLOCKS_PER_SEC);
   printf(":%02d:%02d     %dt\n", s / 60, s % 60, (int)time);
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
   chromlength = numBlocks;
   generation = malloc(sizeof(int *) * popSize);
   for (int i = 0; i < popSize; i++) {
      generation[i] = malloc(sizeof(int) * chromlength);
   }
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
   int input = 0;
   int a = 1;
   while (a != argc) {
      if (!strcmp(argv[a], "-pg")) {
         options->printGen = true;
         options->defaultOptions = false;
      } else if (!strcmp(argv[a], "-t")) {
         options->printTime = true;
         options->defaultOptions = false;
      } else if (!strcmp(argv[a], "-g")) {
         options->greedy = true;
      } else if (!strcmp(argv[a], "-p")) {
         popSize = atoi(argv[a + 1]);
         a++;
      } else if (!strcmp(argv[a], "rand")) {
         if (!input) {
            // instead of determining the output yourself, let every block have
            // a height of min-max
            char path[20];
            snprintf(path, sizeof(path), "rand/%d.in", (int)time(NULL));
            FILE *file = fopen(path, "w");
            if (file == NULL) {
               printf("failure to create %s\n", path);
               exit(EXIT_FAILURE);
            } else {
               printf("created file %s\n", path);
            }
            subsets = atoi(argv[a + 1]);
            numBlocks = atoi(argv[a + 2]);
            blocks = malloc(sizeof(int) * numBlocks);
            fprintf(file, "%d\n", subsets);
            input = 1;
            int randMin = atoi(argv[a + 3]), randMax = atoi(argv[a + 4]) + 1;
            for (int i = 0; i < numBlocks; i++) {
               blocks[i] = (rand() % (randMax - randMin)) + randMin;
               fprintf(file, "%d ", blocks[i]);
            }
            fclose(file);
         }
         a = a + 4;
      } else if (!strcmp(argv[a], "-d")) {
         maxDrift = atoi(argv[a + 1]);
         a++;
      } else {
         if (!input) {
            input = 1;
            FILE *input = fopen(argv[a], "r");
            if (input == NULL) {
               printf("failure to read file\n");
               exit(EXIT_FAILURE);
            }
            fscanf(input, "%d ", &subsets);
            fscanf(input, "%d ", &numBlocks);
            blocks = malloc(sizeof(int) * numBlocks);
            for (int i = 0; i < numBlocks; i++) {
               fscanf(input, "%d ", &blocks[i]);
            }
            fclose(input);
         }
      }
      a++;
   }
   if (!input) {
      (void)fscanf(stdin, "%d ", &subsets);
      (void)fscanf(stdin, "%d ", &numBlocks);
      blocks = malloc(sizeof(int) * numBlocks);
      for (int i = 0; i < numBlocks; i++) {
         (void)fscanf(stdin, "%d ", &blocks[i]);
      }
   }
}

// fill a chromosome with random values
void introduce(int *chromosome) {
   for (int g = 0; g < chromlength; g++) {
      chromosome[g] = towerRand();
   }
}

// mutates 5% of genes by default, changed yes by factor
void mutate(int *chromosome, double factor) {
   int idx = rand() % chromlength;
   int rndm;
   // mutate 5% of the genes, multiplied by factor, and at least 1
   int mutations = (int)((0.05 + factor) * (double)chromlength) + 1;
   for (int i = 0; i < mutations; i++) {
      do {
         rndm = towerRand();
      } while (rndm == chromosome[idx]);
      chromosome[idx] = rndm;
   }
}

void printChrom(int *ranks) {
   printf("Chromosome: ");
   for (int i = 0; i < chromlength; i++) {
      printf("%c", toBase64(generation[ranks[0]][i]));
   }
   printf("\n");
}

void printFullSum(int *ranks, double deviation) {
   int first, sum;
   int max = totalHeight/subsets, min = max;
   for (int s = 0; s < subsets; s++) {
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
      if (sum < min) {
         min = sum;
      }
      if (sum > max) {
         max = sum;
      }
      printf("=%d; dev: %.3f\n", sum,
             (double)sum - ((double)totalHeight / (double)subsets));
   }
   printf("Optimal: %.3f; Average deviation: %.3f; Min: %d; Max: %d\n",
          (double)totalHeight / (double)subsets, deviation, min, max);
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
