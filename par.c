#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "tree.h"

/* par.c contains functions that will be used for the parallelized decision
 * tree construction. We use MPI and design it so that each processor holds
 * the data for one feature plus the labels and a pointer to the sample which
 * can be used as a key. 
 * 
 * Each processor sorts its own data once in the beginning, 
 *
 */

double ParallelSplit(double **data, int n, int first) {

    int rank, p;
    int tag = 21; //Timmy
    int i, j;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int num_features = (D-1)/p;
    if (rank < D % p)
        num_features += 1;

    //////COPYDATA/////////////////
    ALLDATA = MNIST();

    //Allocate pointers to pods 
    Pod **data = malloc(n*sizeof(Pod*));
    //Each process stores array of pods containing feature
    for (i = 0; i < n; ++i)
        data[i] = malloc(sizeof(Pod));

    for (i = 0; i < n; ++i) {
        data[i]->key = i; //key
        data[i]->val = ALLDATA[i][rank]; //feature
        data[i]->label = ALLDATA[i][D-1]; //label
        data[i]->weight = 1./n; //weight
    }
    ///////////////////////////////

    PodSort(data, 0, n-1, 1);
 
    //Get initial counts for positive/negative labels
    int pos = 0;
    double pos_w = 0;//positive weight
    double tot = 0;//total weight
    for (i = first; i < first+n; ++i) {
        tot += data[i]->weight;

        if (data[i]->label > 0){
            pos += 1;
            pos_w += data[i]->weight;
        }
    }
    int neg = n - pos;
    double neg_w = tot - pos_w;

    int row; //best row to split at for particular column/feature
    double threshold; //best threshold to split at for column/feature
    double impurity; //impurity for best split in feature/column

    //get_timestamp(&split_start);
    row = PodWBS(data, n, first, pos_w, tot, &impurity);
    //get_timestamp(&split_stop);
    threshold = data[row]->val;
    //split_time += timestamp_diff_in_seconds(split_start, split_stop);

    //AllReduce to find min impurity and corresponding process (MPI_MINLOC)
    struct {
        double P;
        int R;
    } in, out;

    in.P = impurity;
    in.R = rank;
    
    MPI_Allreduce(in, out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    //If process has min, construct and broadcast list telling which node each point goes to
    int *keys = malloc(n*sizeof(int));
    char *dest = malloc(n*sizeof(char));
    int split_after;
    if (rank == out.R) {
        split_after = row;
        for (i = first; i < row+1; ++i) {
            keys[i] = data[i]->key;
            dest[i] = 0;
        }
        for (i = row+1; i < first+n; ++i) {
            keys[i] = data[i]->key;
            dest[i] = 1;
        }
    }

    MPI_Bcast(keys, n, MPI_INT, out.R, MPI_COMM_WORLD);
    MPI_Bcast(dest, n, MPI_CHAR, out.R, MPI_COMM_WORLD);
    MPI_Bcast(&split_after, 1, MPI_INT, out.R, MPI_COMM_WORLD);
    //Sort pod pointer list into ordered right node and left node
    Pod **newdata = malloc(n*sizeof(Pod*));
    for (i = first; i < first+n; ++i) {
        for (j = 0; j < n; ++j) {
            if (data[i]->index == keys[j]) {
                
        
    
    //Create tree nodes according to split
    
    //Make MPI Barrier, then begin next round; check level, entropy, or purity to decide
    
    

}






