/*
 * reference function to the genetic algorithm version
 */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char **argv) {
    srand(time(NULL));
    for (int i = 0; i < 1000; i++) {
      printf("%d ", rand()%500);
    }
}