// rather simple implementation of the greedy heuristic to solving the number
// partitioning problem

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

int compare(const void *a, const void *b) { return (*(int *)b - *(int *)a); }

int main() {
  int subsetc, amount;
  scanf("%d", &subsetc);
  scanf("%d", &amount);
  
  int *weights = malloc(sizeof(int) * amount);
  int *sums = calloc(sizeof(int), subsetc);

  for (int i = 0; i < amount; i++) {
    scanf("%d ", &weights[i]);
  }
  qsort(weights, amount, sizeof(int), compare);

  int smallestIdx;
  for (int i = 0; i < amount; i++) {
    smallestIdx = 0;
    for (int j = 1; j < subsetc; j++) {
      if (sums[smallestIdx] > sums[j]) {
        smallestIdx = j;
      }
    }
    sums[smallestIdx] += weights[i];
  }
  
  for (int i = 0; i < subsetc; i++) {
    printf("%d: %d\n",i,sums[i]);
  }
  free(sums);
  free(weights);
  return 0;
}
