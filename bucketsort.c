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
#include <math.h>
#include <limts.h>
#include <string.h>
#include <time.h>

typedef struct {
	int count;
	long *bucket_a;
	int bound;
} Bucket;

// Function Prototypes
int valid_sort(const long *array, const int len); 
void merge(long *array, int n, int m);
void serialMergeSort(long *array, int len);
void gen_random_array(long *array, const int len); 
void analyzeSort(long *array, int num_elements, double time, char *type);
long get_random_index(const long *array, int len);
long min(long first, long second);
void p0_setup(long *array_serial, long *array_parallel, int n, int comm_sz,
int *pivots);

int main(int argc, char *argv[]){
	
	int comm_sz; /* Number of processes */
	int my_rank;   /* My process rank*/

	// MPI initializations
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long *array_serial;
    long *array_parallel;
	long *local_array;
	int n = 0;
	int local_n = n/comm_sz;
	int *pivots;

	// Process 0 gets arg and creates arrays with random vals
	//TODO: Error check
	if(my_rank == 0){
		printf("Enter number of elements to sort\n");
		scanf("%d", &n);
	}
	// Broadcast N
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	// Arrays
	array_parallel = malloc(n*sizeof(long));
	local_array = malloc(sizeof(long) * local_n);
	pivots = malloc(sizeof(int) * (comm_sz-1));

	// Process 0 sorts serial array
	if(my_rank == 0) {
		array_serial = malloc(n*sizeof(long));
		p0_setup(array_serial, array_parallel, n, comm_sz, pivots);
	}
	
	// Broadcast pivots and scatter parallel array
	MPI_Bcast(pivots, comm_sz-1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(array_parallel, local_n, MPI_LONG, local_array,
	local_n, MPI_LONG, my_rank, MPI_COMM_WORLD);//Each process has local array

    /*
	 *****correctly calcs and bdcsts pivots, and each process has it's local array.******
     * 
     * TODO:Each process determines which bucket each element of it's local array
     		belongs to. Then send that element to process[i]'s bucket.
     * TODO:Once all processes have all the element assigned to their bucket they
     * 		can create a sorted array of elements in their bucket and send their 
	 		sorted elements back to P0.
     */
	 int i, j;
	 // Create array of "buckets"
	 // The we only make pivot values for the first P-1 processes. The last
	 // process should get anything that is greater than the last pivot.
	 // To get everything, we make the last bucket's bound effectively
	 // infinite.
	 Bucket sim_buckets[comm_sz];
	 for(i = 0; i < comm_sz-1; i++) {
		 sim_buckets[i].bound = pivots[i];
		 sim_buckets[i].a = malloc(sizeof(local_n));
	 }
	 sim_buckets[comm_sz-1].bound = INT_MAX;
	 sim_buckets[comm_sz-1].a = malloc(sizeof(local_n));

	 for(i = 0; i < local_n; i++) { //loop through local elements
		 for(j = 0; j < comm_sz; j++) { // loops through buckets
			 if(local_array[i] <= sim_buckets[j].bound) {
				 Bucket sb = sim_buckets[j];
				 sb.a[sb.count++] = local_array[i];
				 // Break out of bucket loop, go to next element. 
				 break;
			 }
		 }
	 }
	 // TODO Communicate with all the other processes in an offset fashion.
	 // We wrote it down on a sheet of paper. GET THE PAPER. 


    if(my_rank == 0) {
       free(array_serial);
       free(array_parallel);
    }

	MPI_Finalize();
	return 0;
}

void p0_setup(long *array_serial, long *array_parallel, int n, int comm_sz, int *pivots) {

	struct timeval tv1, tv2; //For timing
	// Fill arrays with random values
	gen_random_array(array_serial, n);
	gen_random_array(array_parallel, n);
	
	/* Perform timed serial sort */
	gettimeofday(&tv1,NULL); //Start time
	serialMergeSort(array_serial, n);
	gettimeofday(&tv2,NULL); //End time
	double serial_time = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 +
		(double) (tv2.tv_sec - tv1.tv_sec); 
	analyzeSort(array_serial, n, serial_time, "Serial"); 

	/*** COMPUTE PIVOTS ***/
	// Calc S and create sample array and fill it with random samples
	long num_samples = 10 * comm_sz * (log(n)/log(2));// S
	num_samples = min(n,num_samples);
	long sample_indices[num_samples]; 
	int i;
	srand(time(NULL));
	for(i = 0; i < num_samples; i++) {
		sample_indices[i] = rand() % n;
	}
	// Sort samples
	serialMergeSort(sample_indices, num_samples);
		
	// Find the pivots using pivots[i] = S*(i+1)/P
	for(i = 0; i < comm_sz-1; i++) {
		int w = (num_samples * (i+1)) / comm_sz;
		pivots[i] = sample_indices[w];
		printf("Pivot %d : = %d\n",i,  pivots[i]);
	}
	// pivots now contains the list of pivots
}

/* Finds the minimum of two longs
 */
long min(long first, long second){
	if(first < second)
		return first;
	else
		return second;
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
    long temp[n];
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

void analyzeSort(long *array, int num_elements, double time, char *type){
	if(!valid_sort(array, num_elements)){
		printf("INVALID SORT\n");
	} else{
		printf("======================\n");
		printf("VALID SORT\n");
		//TODO: Add time and complexity stats
		printf("Type: %s\n", type);
		printf("Time: %lf\n", time);
	}
}

