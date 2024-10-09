/*
 * by Rik de Hoop & Nikolaos Trigonis
 * A Genetic Algorithm for the Partition Problem
 * AKA: parprob.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define numBlocks 20
#define chromlength numBlocks
#define popSize 10

int blocks[numBlocks];
bool generation[popSize][chromlength];

int readHeights() {}
int fillChromosome() {}
int mutate() {}
int crossOver() {}
int heightOfTower(bool *chromosome, bool b) {}
int heightDif(bool *chromosome) {}

int main(int argc, char **argv) {
  srand(time(NULL));
  readHeights();


  // for(int i = 0; i < numBlocks; i++) {
  // (void)scanf("%d ", &blocks[i]);
  // }

  // are you allowed to initialize with (a) heuristic(s)?

}

int heightOfTower(bool *chromosome, bool b) {
  int height = 0;
  if (b) {
    for (int i = 0; i < chromlength; i++) {
      if (chromosome[i]) {
        height += blocks[i];
      }
    }
  } else {
    for (int i = 0; i < chromlength; i++) {
      if (!chromosome[i]) {
        height += blocks[i];
      }
    }
  }
  return height;
}

int heightDif(bool *chromosome) {
  int dif = 0;
  for (int i = 0; i < chromlength; i ++) {
    if (chromosome[i]) {
      dif += blocks[i];
    } else {
      dif -= blocks[i];
    }
  }
  return abs(dif);
}
