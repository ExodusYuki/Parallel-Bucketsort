/* Melinda Grad
 * Quentin Fulsher
 * PA3
 */
#include <sys/time.h>
#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <time.h>

// Function Prototypes
int valid_sort(const long *array, const int len); 
void merge(long *array, int n, int m);
void serialMergeSort(long *array, int len);
void print_array(const long *array, const int len);
void gen_random_array(long *array, const int len); 

int main(int argc, char *argv[]){
	
	int comm_size; /* Number of processes */
	int my_rank;   /* My process rank*/
	struct timeval tv1, tv2; //For timing

	// MPI initializations
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Process 0 gets arg and creates arrays with random vals
	if(my_rank == 0){
		if(argc != 2){
			printf("Program needs 1 arg: <number of elements>\n");
			exit(1);
		}	
		int n = strtol(argv[1], NULL, 10);

		// Create two arrays with same random values
		long *array_serial = malloc(n*sizeof(long));
		long *array_parallel = malloc(n*sizeof(long));
		gen_random_array(array_serial, n);
		gen_random_array(array_parallel, n);
		
		/* Perform timed serial sort */
		gettimeofday(&tv1,NULL); //Start time
		serielMergeSort(array_parallel, n);
		gettimeofday(&tv2,NULL); //End time
		double serial_time = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 +
			(double) (tv2.tv_sec - tv1.tv_sec); 
		analyzeSort(array_serial, n, serial_time, "Serial"); 


	}
	//TODO: Print Data (#7 on project writeup)

	MPI_Finalize();
	free(array_serial);
	free(array_parallel);
	return 0;
}

/*
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

/* 
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

/*
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

/**
* Validates the sort was successful.
* returns true if valid, false if invalid
* (Uses stdbool types)
*/
int validSort(const int *array, const int len) {
	int i;
	for (i = 1; i < len; i++) {
		if(array[i-1] > array[i]) {
			return false;
		}
	}
	 return true;
}

void analyzeSort(int *array, int num_elements, double time, char *type){
	if(!validSort(array, num_elements)){
		printf("INVALID SORT\n");
	} else{
		printf("======================\n");
		printf("VALID SORT\n");
		//TODO: Add time and complexity stats
	}
}

