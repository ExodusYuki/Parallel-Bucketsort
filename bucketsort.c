/* Melinda Grad
 * Quentin Fulsher
 * PA3
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

// Function Prototypes
int valid_sort(const int *array, const int len); 
void merge(int *array, int n, int m);
void serialMergeSort(int *array, int len);
void print_array(const int *array, const int len); 
int main(int argc, char *argv[]){

	// Arg is number of elements in the array to be sorted
	// Assume it is a multiple of the number of processes
	if(argc != 2){
		printf("Program needs 1 arg: <number of elements>\n");
		exit(1);
	}
	// Process 0 does this
	int num_elements = strtol(argv[1], NULL, 10);
	//TODO: Process 0 creates 2 arrays (long) with same
	// random values one seriel, one parallel
	

}

/* Added from PA1
 * Regular seriel mergesort implementation
 */
void serialMergeSort(int *array, int len){
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
void merge(int *array, int n, int m){
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
int valid_sort(const int *array, const int len) {
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
void print_array(const int *array, const int len) {
	int i;
	for (i = 0; i < len; i++) {
		printf("%d, ", array[i]);
	}
	printf("\n");
}
