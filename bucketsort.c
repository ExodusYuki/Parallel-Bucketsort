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
#include <limits.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {
	int count;
	long *a;
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
Bucket *createBuckets(int comm_sz, int *pivots, int local_n, long *local_array);
int sendRecvBuckets(int my_rank, int comm_sz, Bucket *sim_buckets, long *recv_buff);
void updateRecvPartner(int *recv_partner, int comm_sz);
void updateSendPartner(int *send_partner, int comm_sz);
void printAllBuckets(Bucket *sim_buckets, int my_rank, int len);
void printAllArray(long *array, int my_rank, int len);
void printArray(long *array, int len);
void printBuckets(Bucket *sim_buckets, int num_buckets);

int main(int argc, char *argv[]){
	
	int comm_sz; /* Number of processes */
	int my_rank;   /* My process rank*/	
	long *array_serial;
    long *array_parallel;
	long *local_array;
	int n;
	int local_n;
	int *pivots;
    struct timeval tv1, tv2;

	// MPI initializations
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Process 0 gets arg and creates arrays with random vals
	//TODO: Error check
	if(my_rank == 0){
		printf("Enter number of elements to sort\n");
		scanf("%d", &n);
	}
	
	// Broadcast N
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	local_n = n/comm_sz;
	
	// Allocate arrays
	array_parallel = malloc(sizeof(long) * n);
	local_array = malloc(sizeof(long) * local_n);
	pivots = malloc(sizeof(int) * (comm_sz-1));
		
	if(my_rank == 0) {
		array_serial = malloc(n*sizeof(long));
		p0_setup(array_serial, array_parallel, n, comm_sz, pivots);
        gettimeofday(&tv1, NULL);
	}

    // Broadcast pivots
    MPI_Bcast(pivots, comm_sz-1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(array_parallel, // Distribute the array
        local_n,					// Number each proceess handles
        MPI_LONG, 					// TYPE
        local_array,				// Receive Buffer
        local_n,					// Receive count
        MPI_LONG, 					//Type
        0,							//Root
        MPI_COMM_WORLD
    );
    // Each now has local array

	// Create array of "buckets"
	Bucket *sim_buckets = createBuckets(comm_sz, pivots, local_n, local_array);

	long recv_buff[n];
 	int recv_buff_sz = sendRecvBuckets(my_rank, comm_sz, sim_buckets, recv_buff);
  	serialMergeSort(recv_buff, recv_buff_sz);

	/* Process 0 gathers the buckets-- see book pg. 113*/
    long *final_array = NULL;
	if(my_rank == 0){
		final_array = calloc(sizeof(long), n);
	}

    /*
     * Put bucket stuff in array, then recv from other processes into the final array
     */
    if(my_rank == 0) {
        int off = 0;
        for(off = 0; off < recv_buff_sz; off++) {
            final_array[off] = recv_buff[off];
        }
        int i;
        int incoming_size = 0;
        for(i = 1; i < comm_sz; i++) {
            MPI_Recv(&incoming_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            MPI_Recv(final_array+off, incoming_size, MPI_LONG, i, 0, MPI_COMM_WORLD, NULL);
            off += incoming_size;
        }
    } else {
        MPI_Send(&recv_buff_sz, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(recv_buff, recv_buff_sz, MPI_LONG, 0, 0, MPI_COMM_WORLD);
    }
		
    if(my_rank == 0) {
        gettimeofday(&tv2, NULL);
        double parallel_time = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 +
            (double) (tv2.tv_sec - tv1.tv_sec); 
        analyzeSort(final_array, n, parallel_time, "Parallel");
    }

    free(sim_buckets);
    if(my_rank == 0) {
       free(array_serial);
       free(array_parallel);
    }

	MPI_Finalize();
	return 0;
}


/*TODO: Check and finish me
 * Communicate with all the other processes in an offset fashion.
 * We wrote it down on a sheet of paper. GET THE PAPER. 
 */
int sendRecvBuckets(int my_rank, 
        int comm_sz, Bucket *sim_buckets, long *recv_buff){

	int i = 0;
    int recv_size = 0;
    int num_recvd = 0;
	int recv_partner = my_rank;
	int send_partner = my_rank;
    updateRecvPartner(&recv_partner, comm_sz);
    updateSendPartner(&send_partner, comm_sz);
	for(i = 0; i < comm_sz; i++){
        // Send size of the bucket to send_partner,
        // recv size of incoming buffer from recv_partner
		MPI_Sendrecv(
            &sim_buckets[send_partner].count, // Send size of outgoing array
            1,                              // We are sending one value
            MPI_INT,                        // Type = int
            send_partner,                   // Sending it to our send partner
            0,                              // TAG
            &recv_size,                     // We want to recv into the recv_size variable
            1,                              // we get one value
            MPI_INT,                        // Also of type int
            recv_partner,                   // From our recv partner
            0,                              // TAG
            MPI_COMM_WORLD,                 // Communicator
            MPI_STATUS_IGNORE               // Status
        );

		MPI_Sendrecv(
            sim_buckets[send_partner].a,     // Send the array
            sim_buckets[send_partner].count,// We want to send the whole array
            MPI_LONG,                       // TYPE = LONG
            send_partner,                   // Send to send partner
            0,                              // TAG
            recv_buff+num_recvd,            // Recv into recv_buff+number of elements already recvd
            recv_size,                      // Recv amount of elements = to recv_size
            MPI_LONG,                       // TYPE = LONG
            recv_partner,                   // recv from recv_partner
            0,                              // TAG
            MPI_COMM_WORLD,                 // Communicator
            MPI_STATUS_IGNORE               // Status
        );

        num_recvd += recv_size;
        updateRecvPartner(&recv_partner, comm_sz);
        updateSendPartner(&send_partner, comm_sz);
	}
    return num_recvd;
}

void updateRecvPartner(int *recv_partner, int comm_sz) {
    *recv_partner -= 1;
    if(*recv_partner <  0) {
        *recv_partner = (comm_sz - 1);
    }
}

void updateSendPartner(int *send_partner, int comm_sz) {
    *send_partner = (*send_partner + 1) % comm_sz;
}
    
/*
 Create array of "buckets"
 The we only make pivot values for the first P-1 processes. The last
 process should get anything that is greater than the last pivot.
 To get everything, we make the last bucket's bound effectively infinite.
*/
Bucket *createBuckets(int comm_sz, int *pivots, int local_n, long *local_array){
	 int i, j;

	 Bucket *sim_buckets = calloc(sizeof(Bucket), comm_sz);
	 for(i = 0; i < comm_sz-1; i++) {
		 sim_buckets[i].bound = pivots[i];
		 sim_buckets[i].a = calloc(sizeof(long), local_n);
	 }
	 sim_buckets[comm_sz-1].bound = INT_MAX;
	 sim_buckets[comm_sz-1].a = calloc(sizeof(long), local_n);
     
	 // Figure out who goes in each bucket
	 for(i = 0; i < local_n; i++) { //loop through local elements
		 for(j = 0; j < comm_sz; j++) { // loops through buckets
			 if(local_array[i] <= sim_buckets[j].bound) {
			    // This is just to assign each element into their bucket.
                // All elements need to be in a bucket before everyone sends
                int local_count = sim_buckets[j].count;
                sim_buckets[j].a[local_count]= local_array[i];
                sim_buckets[j].count += 1;
                // break from inner loop vvi
                break;
			 }

		 }
	 }
     return sim_buckets; // TODO Free this
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
		sample_indices[i] = array_parallel[rand() % n];
	}
	// Sort samples
	serialMergeSort(sample_indices, num_samples);
		
	// Find the pivots using pivots[i] = S*(i+1)/P
	for(i = 0; i < comm_sz-1; i++) {
		int w = (num_samples * (i+1)) / comm_sz;
		pivots[i] = sample_indices[w];
	}
}

void printAllBuckets(Bucket *sim_buckets, int my_rank, int len) {
   	if(my_rank == 0) {
		printBuckets(sim_buckets, len);
	} else if(my_rank == 1) {
        sleep(1);
		printBuckets(sim_buckets, len);
    } else if(my_rank == 2) {
        sleep(2);
		printBuckets(sim_buckets, len);
    } else if(my_rank == 3) {
        sleep(3);
		printBuckets(sim_buckets, len);
    }
}

void printAllArray(long *array, int my_rank, int len) {
   	if(my_rank == 0) {
        printf("Array for #%d\n", my_rank);
		printArray(array, len);
	} else if(my_rank == 1) {
        sleep(1);
        printf("Array for #%d\n", my_rank);
		printArray(array, len);
    } else if(my_rank == 2) {
        sleep(2);
        printf("Array for #%d\n", my_rank);
		printArray(array, len);
    } else if(my_rank == 3) {
        sleep(3);
        printf("Array for #%d\n", my_rank);
		printArray(array, len);
        printf("====================\n");
    }
}

void printArray(long *array, int len) {
    int i;
    for(i = 0; i < len; i++) {
        printf("%ld ", array[i]);
    }
    printf("\n");
}

void printBuckets(Bucket *sim_buckets, int num_buckets) {
    int i;
    for(i = 0; i < num_buckets; i++) {
        Bucket sb = sim_buckets[i];
        printf("Bucket #%d:\n", i+1);
        printf("\tcount: %d\n", sb.count);
        printf("\tbound: %d\n", sb.bound);
        int j;
        printf("\tarray:");
        for(j = 0; j < sb.count; j++) {
           printf("%ld ", sb.a[j]);
        }
        printf("\n");
        printf("=================\n");
    }
    printf("\n");
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
            printf("%ld > %ld\n", array[i-1], array[i]);
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

