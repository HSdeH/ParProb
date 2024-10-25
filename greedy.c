// rather simple implementation of the greedy heuristic to solving the number
// partitioning problem

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define AMOUNT 100
#define K 3

int compare(const void *a, const void *b) { return (*(int *)b - *(int *)a); }

int main() {
  int *weights = malloc(sizeof(int) * AMOUNT);
  int *sums = calloc(sizeof(int), K);

  for (int i = 0; i < AMOUNT; i++) {
    scanf("%d ", &weights[i]);
  }
  qsort(weights, AMOUNT, sizeof(int), compare);

  int smallestIdx;
  for (int i = 0; i < AMOUNT; i++) {
    smallestIdx = 0;
    for (int j = 1; j < K; j++) {
      if (sums[smallestIdx] > sums[j]) {
        smallestIdx = j;
      }
    }
    sums[smallestIdx] += weights[i];
  }
  
  for (int i = 0; i < K; i++) {
    printf("%d: %d\n",i,sums[i]);
  }
  free(sums);
  free(weights);
  return 0;
}