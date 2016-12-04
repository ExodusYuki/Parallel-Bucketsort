/* Melinda Grad
 * Quentin Fulsher
 * PA3
 */
#include <sys/time.h>
#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
//#include <mpi.h>
#include <string.h>
#include <time.h>

// Function Prototypes
int valid_sort(const long *array, const int len); 
void merge(long *array, int n, int m);
void serialMergeSort(long *array, int len);
void print_array(const long *array, const int len);
void gen_random_array(long *array, const int len); 

int main(int argc, char *argv[]){
	//TODO: Add timing
	struct timeval tv1, tv2; //For timing

	if(argc != 2){
		printf("Program needs 1 arg: <number of elements>\n");
		exit(1);
	}
	//TODO: Process 0 gets arg and creates arrays with random vals
	// Number of elements in array
	int n = strtol(argv[1], NULL, 10);

	// Create two arrays with same random values
	long *array_serial = malloc(n*sizeof(long));
	long *array_parallel = malloc(n*sizeof(long));
	gen_random_array(array_serial, n);
	gen_random_array(array_parallel, n);

	// Debug
	print_array(array_serial, n);
	print_array(array_parallel, n);

	return 0;
}

/*Added from PA1
 *Takes an array and generates random numbers up to len
 */
void gen_random_array(long *array, const int len) {
	int i;
	//Make sure our random values aren't the same everytime
	srand(time(NULL));
	//Fill the array with random ints
	for (i = 0; i < len; i++) {
		array[i] = (rand() % 100);
	}
}

/* Added from PA1
 * Regular seriel mergesort implementation
 */
void serialMergeSort(long *array, int len){
	// Base Case
	if(len < 2){
		return;
	}
	int mid = len/2;
	//Mergesort the left half
	serialMergeSort(array, mid);
	//Mergesort the right half
	serialMergeSort(array+mid, len-mid);
	merge(array, len, mid);
}

/*TODO: Change to work with new arrays
 * Added from PA1
 * Regular merge func to be used with seriel mergesort
 */
void merge(long *array, int n, int m){
	int i, j, k;
	for(i = 0, j = m, k = 0; k < n; k++){
		if(j ==n){
			//If right is out of elements
			temp[k] = array[i++];
		} else if(i == m){
			// If left is out of elements
			temp[k] = array[j++];
		} else if(array[j] < array[i]){
			temp[k] = array[j++];
		} else{
			temp[k] = array[i++];
		}
	}
	// Copy numbers back to parent array
	for(i = 0; i<n; i++){
		array[i] = temp[i];
	}
}

/**
 * Added from PA1
 * Validates the sort was successful,
 * returns true if valid, false if invalid
 * (Uses std bool types)
 */
int valid_sort(const long *array, const int len) {
	int i;
	for(i = 1; i< len; i++){
		if(array[i-1] > array[i]){
			return false;
		}
	}
	return true;
}

/**
 *Simple function that prints an int array of length len on one line.
 * */
void print_array(const long *array, const int len) {
	int i;
	for (i = 0; i < len; i++) {
		printf("%ld, ", array[i]);
	}
	printf("\n");
}
