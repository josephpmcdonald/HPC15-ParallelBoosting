#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Notes: Use main to test functions.
 * 
 * TODO:
 * In BestSplit, take care of case when i=n before computing impurity, else division by 0
 *
 *
 */


#define GINI(pos, n) (2.*(pos)/(n)*((n)-(pos))/(n))
#define ENTROPY(pos, n) ((pos)/(n)*log((pos)/(n)) + ((n)-(pos))/(n)*log(((n)-(pos))/(n)))

const int classes = 2;
const int N = 1000;
const int d = 10;

struct Node {
    struct Node *parent;
    struct Node *right;
    struct Node *left;
    int index;
    double threshold;
    double class;
};


void Sort(float data[][d], int first, int last, int a) {

/* Sort and Partition implements quicksort and sorts a 2-dimensional
 * array on a particular index a between indices first and last. Sort places
 * entries matching the pivot next to the pivot in the sorted array to minimize
 * recursive Sort calls in the case of repeated values.
 *
 * data  = array containing sample data and labels to be sorted
 * first = first unsorted row that Sort may move
 * last  = last unsorted row that Sort may move
 * a     = index to sort on
 *
 */

    if (first >= last)
        return;

    //Save the ends of the unsorted portions in vector.
    int ends[2];

    Partition(data, first, last, a, ends);
    Sort(data, first, ends[0], a);
    Sort(data, ends[1], last, a);
    
    return;
}


void Partition(float data[][d], int first, int last, int a, int ends[]) {

/* Partition implements the partition function in quicksort. It groups
 * rows that match the pivot value together to minimize recursive calls
 * to Sort. 
 *
 * data  = array containing sample data and labels to be sorted
 * first = first unsorted row that Sort may move
 * last  = last unsorted row that Sort may move
 * a     = index to sort on
 * ends  = vector with last index of lower unsorted block and first 
 *         index of higher unsorted block
 *
 */

    double pivot = data[last][a];
    int i = first;
    int j = last-1;
    int k = last;
    int ind;
    double holder;
    while (i < j) {
        if (data[i][a] < pivot)
            ++i;
        else if (data[i][a] > pivot) {
            //swap i with j, reduce j
            for (ind = 0; ind < d; ++ind) {
                holder = data[j][ind];
                data[j][ind] = data[i][ind];
                data[i][ind] = holder;
            }
            --j;
        }
        else { //(data[i][a] == pivot)
            //i goes to k-1, k-1 to j, and j to i, reduce j, reduce k
            for (ind = 0; ind < d; ++ind) {
                holder = data[i][ind];
                data[i][ind] = data[j][ind];
                data[j][ind] = data[k-1][ind];
                data[k-1][ind] = holder;
            }
            --j;
            --k;
        }
    }

    //Now i = j.
    if (data[i][a] < pivot)
        ++i;
    else if (data[i][a] > pivot)
        --j;
    else { // (data[i][a] == pivot) 
        for (ind = 0; ind < d; ++ind) {
            holder = data[i][ind];
            data[i][ind] = data[k-1][ind];
            data[k-1][ind] = holder;
        }
        --j;
        --k;
    }

    //Now i = j+1. j and below are less than pivot, i and above are greater.
    //Last-k+1 entries match pivot. Swap to the dividing region, starting at i.
    int matches = last-k+1;
    for (l = 0; l < matches; ++l) {
        for (ind = 0; ind < d; ++ind) {
            holder = data[i+l][ind];
            data[i+l][ind] = data[k+l][ind];
            data[k+l][ind] = holder;
        }
    }
    
    ends[0] = i-1;
    ends[1] = i + matches;

    return;
}


int BestSplit(float data[][d], int n, int first, int col, int pos, double *impurity) {

/* Returns the row/index of the table with the least impurity after splitting
 * for fixed column/feature col. Partition rows up to and including that index
 * from everything afterwards.
 *
 * data = array of data sorted on index a 
 * n    = length of table (# of rows/samples)
 * col    = sorting/splitting feature/column of data
 * pos  = number of positive labels
 *
 * Thus:
 * (pos - lpos) = right pos count
 * i = total left count
 * n - i = total right count
 * 
 *
 */

    int neg = n - pos;
    int lpos = 0;

    int argmin = n - 1;//start with the whole node
    double threshold;
    double threshmin;
    double P;
    double Pmin = GINI(pos, n);//initial impurity of node

    //Tabulate impurity for each possible threshold split    
    int i = 0;
    while (i < n) {

        threshold = data[first+i][col];

        while (data[first+i][col] == threshold) {
            if (data[first+i][d-1] == 1)
                lpos += 1;

            ++i;
        }

        //Note that points on left = i, right = n-i

        /*
        If i=n, this is the whole node and the impurity is the initial
        which is already done. i=n would cause error below.
        */
        if (i < n) {
            P = GINI(lpos, i)*i/n + GINI(pos-lpos, n-i)*(n-i)/n;

            //Save threshold/index with min impurity
            if (P < Pmin) {
                Pmin = P;
                argmin = i - 1;
                threshmin = threshold; 
            }
        }
    }

    //Save the minimum impurity to compare against other indices
    *impurity = Pmin;
    return argmin;
}


void SplitNode(struct Node *node, float data[][d], int n, int first, int level) {

/* Creates two branches of the decision tree on the array data. End condition
 * creates leaf if the purity of the node is small or if there are few
 * samples on the branch of node
 *
 * node  = pointer to node in decision tree
 * data  = table of unsorted data with features and labels (with last
 *         column as the label (data[i][d-1]))
 * n     = length of table (# of rows/samples) on branch of node
 * first = first index of samples on branch of node
 */

    //Get initial counts for positive/negative labels
    int i;
    int pos;
    for (i = 0; i < n; ++i) {
        if (data[first+i][d-1] == 1.)
            pos += 1;
    }
    int neg = n - pos;

    //Declare class for node in case of pruning on child
    if (pos > neg)
        node->class = 1;
    else if (pos < neg)
        node->class = -1;
    else if (node->parent)
        node->class = node->parent->class;
    else {
        fprint("Root node is evenly balanced.\n");
        node->class = 0;
    }

    //If branch is small or almost pure, make leaf
    if (n < 6 || GINI(pos, n) < 0.01 || level == 6) {
        node->left = NULL;
        node->right = NULL;
        return;
    }


    int col;
    int row; //best row to split at for particular column/feature
    int localrow; //first + localrow = row; receives BestSplit which returns integer in [-1, n-1]
    double threshold; //best threshold to split at for column/feature
    double impurityaddress;
    double *impurity = &impurityaddress; //pointer to impurity for feature
    int bestcol; //feature with best split
    int bestrow; //best row to split for best feature
    double bestthresh; //threshold split for best feature (data[bestrow][bestcol])
    double Pmin = 100; //minimum impurity seen so far

    //Sort table. Then find best column/feature, threshold, and impurity
    for (col = 0; col < d-1; ++col) {
        Sort(data, first, first+n-1, col);
        localrow = BestSplit(data, n, first, col, pos, impurity);
        row = first + localrow;
        threshold = data[row][col];

        //If current column has better impurity, save col, thresh, and Pmin
        if (*impurity < Pmin) {
            bestcol = col;
            bestrow = row;
            bestthresh = threshold;
            Pmin = *impurity;
        }
    }

    //If splitting doesn't improve purity (best split is at the end) stop
    if (bestrow == first+n-1) {
        node->left = NULL;
        node->right = NULL;
        return;
    }

    Sort(data, first, first+n-1, bestcol);

    //Quick possibly unnecessary check just to make sure code works correctly
    if (bestthresh != data[bestrow][bestcol])
        fprint("ERROR IN SplitNode: TREE BUILDING MESSED UP.\n");

    //For feature, threshold with best impurity, save to node attributes
    node->index = bestcol;
    node->threshold = bestthresh;

    //Create right and left children
    struct Node right;
    struct Node left;
    right.parent = node;
    left.parent = node;
    node->right = &right;
    node->left = &left;

    int first_r = bestrow+1;
    int n_l = first_r - first;
    int n_r = n - n_l;

    SplitNode(left, data, n_l, first, level+1);
    SplitNode(right, data, n_r, first_r, level+1);

    return;
}



/*
    double *ldata = malloc(argmin * sizeof(double*));
    double *rdata = malloc((n - argmin) * sizeof(double*));
    double *lmatrix = malloc(argmin * d * sizeof(double));
    double rmatrix = malloc((n - argmin) * d * sizeof(double));
    for (i = 0; i < argmin; ++i)
        ldata[i] = &lmatrix[i*d];
    for (i = 0; i < n-argmin; ++i)
        rdata[i] = &rmatrix[i*d];

    int j = 0;
    for (i = 0; i < argmin; ++i) {
        for (j = 0; j < d; ++j)
            ldata[i][j] = data[i][j];
    }

    for (i = 0; i < n - argmin; ++i) {
        for (j = 0; j < d; ++j)
            rdata[i][j] = data[i + argmin][j];
    }
*/




struct node* CountSplit(float data[][d], int n, int a, node *) {

/* data = node data
 * n = length of data
 * a = index to split on
 *
 * Table computes/stores impurity for splitting with data <= threshold
 * counts[i][0] = ith threshold
 * counts[i][1] = left pos countcounts[i][2] = left neg count = i+1 - left pos count
 * counts[i][3] = impurity (Gini: 2*p*(1-p), entropy: p*log p + (1-p)log(1-p))
 *
 * Thus:
 * (pos - counts[i][1]) = right pos count
 * (neg - counts[i][2]) = right neg count
 * i+1 = total left count
 * n-(i+1) = total right count
 *
 */
    
    //
    //Counts stores pointer to each row of matrix
    double *matrix = calloc(n*3, sizeof(double));
    if (matrix)
        double *counts[] = malloc(n*sizeof(double*));

    int pos = 0;
    int neg = 0;
    int lpos = 0;
    int left = 0;
    int right = n;
    double lp;
    double rp;
    double purity;
    double argmin;
    double threshmin;

    //Initialize pointers in counts
    int i = 0;
    for (i = 0; i < n; ++i)
        counts[i] = &matrix[3*i];

    //Get the initial counts for positive/negatives.
    for (i = 0; i < n; ++i) {
        if (data[d*i + d-1] == 1.)
            pos += 1;
    }

    neg = n - pos;
    
    i = 0;
    int ind = 0;
    double threshold;
    double Pmin = 100.; //Upper bound on purity
    int argmin;
    while (i < n) {

        threshold = data[i][a];
        counts[ind][0] = threshold;

        while (data[i][a] == threshold) {
            if (data[i][d-1] == 1) {
                counts[ind][1] += 1;
                lpos += 1;
            }

            ++i;
            //Now left = i

            //else lneg = i - lpos

//            else if (data[i][d-1] == -1)
//                counts[ind][2] += 1;
//            else
//                printf("Data not clean.\n");

        }

        //Note that points on left = i, right = n-i
        lp = counts[ind][1]/i;
        //lp = (double) lpos/i;
        rp = (double) (pos - counts[ind][1])/(n - i);
        //rp = (double) (pos - lpos)/(n - i);

        //P = (left/n)*2*lp*(1 - lp) + (right/n)*2*rp*(1 - rp) and counts[ind][1]/n = (left/n)*lp
        counts[ind][2] = 2*(counts[ind][1]/n)*(1 - lp) + 2*(n - i)/n*rp*(1 - rp);
        purity = 2.*lpos/n*(1 - lp) + 2.*(pos - lpos)/n*(1 - rp);

        if (purity < Pmin) {
            Pmin = purity;
            argmin = i; //Note this is the index of the first right point
            threshmin = threshold;
        }

        //++ind;
    }
    
    for (i = 0; i < argmin; ++i) {
        

    }


    return 0;
}






int main(int argc, char *argv[]) { 

    return 0;
}






