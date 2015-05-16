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

 
double ParallelSplit(Node *node, double **data, int n, int first, int level, int rank) {

    int max_level = 2;
    int min_points = 6;

    node->left = NULL;
    node->right = NULL;
    node->index = -1;

    //int tag = 21; //Timmy
    int i, j;
    
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

    PodSort(data, 0, n-1, 1);

    int row; //best row to split at for particular column/feature
    double threshold; //best threshold to split at for column/feature
    double impurity; //impurity for best split in feature/column

    //get_timestamp(&split_start);
    row = PodWBS(data, n, first, pos_w, tot, &impurity);
    //get_timestamp(&split_stop);
    threshold = data[row]->val;
    //split_time += timestamp_diff_in_seconds(split_start, split_stop);

    //Save impurity and process rank to structure
    struct {
        double P;
        int R;
    } in, out;

    in.P = impurity;
    in.R = rank;
    
/* AllReduce to find min impurity and corresponding process (MPI_MINLOC), then
 * receive best row and hence size of next left/right nodes
 */
    MPI_Allreduce(in, out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    MPI_Bcast(&row, 1, MPI_INT, out.R, MPI_COMM_WORLD);
    printf("rank %d row %d\n", rank, row);

    //If splitting doesn't improve purity (best split is at the end) stop
    if (bestrow == first+n-1) {
        //printf("no improvement\n");
        return;
    }

    //For min processor, construct and broadcast list telling which node each point goes to
    int *keys = malloc((row+1-first)*sizeof(int));
    if (rank == out.R) {
        for (i = 0; i < row+1-first; ++i)
            keys[i] = data[first+i]->key;
    }
    MPI_Bcast(keys, n, MPI_INT, out.R, MPI_COMM_WORLD);

    //Sort pod pointer list into ordered right node and left node
    Pod **holder = malloc(n*sizeof(Pod*));
    int l_ind = 0;
    int r_ind = row+1-first;
    for (i = first; i < first+n; ++i) {
        for (j = 0; j < row+1-first; ++j) {
            if (keys[j] == data[i]->index) {
                holder[l_ind] = data[i];
                l_ind++;
                break;
            }
        }
        if (j == row+1-first) {
            holder[r_ind] = data[i];
            r_ind++;
        }
    }

    for (i = 0; i < n; ++i)
        data[first+i] = holder[i];

    free(keys);
    free(holder);

    //Create right and left children
    Node *l = malloc(sizeof(Node));
    Node *r = malloc(sizeof(Node));
    l->parent = node;
    r->parent = node;
    l->right = NULL;
    l->left = NULL;
    r->right = NULL;
    r->left = NULL;
    
    node->left = l;
    node->right = r;
    
    //Make MPI Barrier, then begin next round; check level, entropy, or purity to decide
    
    
}






