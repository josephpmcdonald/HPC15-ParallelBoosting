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

struct DataPod{
    int index;
    double val;
    double label;
    double weight;
}

typedef struct DataPod Pod;

double ParSplit(double **data, int n) {


    int rank, p, tag=0;
    int i;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int my_features = (D-1)/p;
    if (rank < D % p)
        my_features += 1;

    //////COPYDATA/////////////////
    ALLDATA = MNIST();

    //Allocate pointers to pods 
    Pod **data = malloc(n*sizeof(Pod*));
    //Each process stores array of pods containing feature
    for (i = 0; i < n; ++i)
        data[i] = malloc(sizeof(Pod));

    for (i = 0; i < n; ++i) {
        data[i]->index = i; //index
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
    for (i = 0; i < n; ++i) {
        tot += data[first+i]->weight;

        if (data[first+i]->label > 0){
            pos += 1;
            pos_w += data[first+i]->weight;
        }
    }
    int neg = n - pos;
    double neg_w = tot - pos_w;

    int row; //best row to split at for particular column/feature
    double threshold; //best threshold to split at for column/feature
    double impurity; //impurity for best split in feature/column

    //get_timestamp(&split_start);
    row = first + PodWBS(data, n, first, pos_w, tot, &impurity);
    //get_timestamp(&split_stop);
    threshold = data[row][1];
    //split_time += timestamp_diff_in_seconds(split_start, split_stop);

    

}






