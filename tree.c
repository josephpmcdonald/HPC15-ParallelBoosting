#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tree.h"


/* Notes: Use main to test functions.
 * 
 * TODO: Merge
 *
 *
 */


static const int classes = 2;
static const int N = 1000;

void MergeSort(double data[][D], int first, int last, int a) {

    if (first >= last)
        return;

    int mid = (last-first)/2;
    MergeSort(double data[][D], int first, first+mid, int a);
    MergeSort(double data[][D], int first+mid+1, int last, int a);

    Merge(double data[][D], int first, int first+mid, int last, int a);
    return;
}

void Merge(double data[][D], int first, int mid, int last, int a) {
    //TODO

    return;
}

void Sort(double data[][D], int first, int last, int a) {

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


void Partition(double data[][D], int first, int last, int a, int ends[]) {

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
    int l;
    int ind;
    double holder;
    while (i < j) {
        if (data[i][a] < pivot)
            ++i;
        else if (data[i][a] > pivot) {
            //swap i with j, reduce j
            for (ind = 0; ind < D; ++ind) {
                holder = data[j][ind];
                data[j][ind] = data[i][ind];
                data[i][ind] = holder;
            }
            --j;
        }
        //(data[i][a] == pivot)
        else { 
            //i goes to k-1, k-1 to j, and j to i, reduce j, reduce k
            for (ind = 0; ind < D; ++ind) {
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
    //(data[i][a] == pivot)
    else {
        for (ind = 0; ind < D; ++ind) {
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
        for (ind = 0; ind < D; ++ind) {
            holder = data[i+l][ind];
            data[i+l][ind] = data[k+l][ind];
            data[k+l][ind] = holder;
        }
    }
    
    ends[0] = i-1;
    ends[1] = i + matches;

    return;
}


int BestSplit(double data[][D], int n, int first, int col, int pos, double *impurity) {

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
            if (data[first+i][D-1] == 1)
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


void SplitNode(struct Node *node, double data[][D], int n, int first, int level) {

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

    node->left = NULL;
    node->right = NULL;

    //Get initial counts for positive/negative labels
    int i;
    int pos;
    for (i = 0; i < n; ++i) {
        if (data[first+i][D-1] == 1.)
            pos += 1;
    }
    int neg = n - pos;

    //Declare class for node in case of pruning on child
    if (pos > neg)
        node->label = 1;
    else if (pos < neg)
        node->label = -1;
    else if (node->parent)
        node->label = node->parent->label;
    else {
        printf("Root node is evenly balanced.\n");
        node->label = 0;
    }


    //If branch is small or almost pure, make leaf
    if (n < 6 || GINI(pos, n) < 0.01 || level == 6)
        return;

    int col;
    int row; //best row to split at for particular column/feature
    int localrow; //first + localrow = row; receives BestSplit which returns integer in [-1, n-1]
    double threshold; //best threshold to split at for column/feature
    double impurity; //impurity for best split in feature/column
    int bestcol; //feature with best split
    int bestrow; //best row to split for best feature
    double bestthresh; //threshold split for best feature (data[bestrow][bestcol])
    double Pmin = 100; //minimum impurity seen so far

    //Sort table. Then find best column/feature, threshold, and impurity
    for (col = 0; col < D-1; ++col) {
        Sort(data, first, first+n-1, col);
        localrow = BestSplit(data, n, first, col, pos, &impurity);
        row = first + localrow;
        threshold = data[row][col];

        //If current column has better impurity, save col, thresh, and Pmin
        if (impurity < Pmin) {
            bestcol = col;
            bestrow = row;
            bestthresh = threshold;
            Pmin = impurity;
        }
    }

    //If splitting doesn't improve purity (best split is at the end) stop
    if (bestrow == first+n-1)
        return;


    Sort(data, first, first+n-1, bestcol);

    //Quick possibly unnecessary check just to make sure code works correctly
    if (bestthresh != data[bestrow][bestcol])
        printf("ERROR IN SplitNode: TREE BUILDING MESSED UP.\n");

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

    SplitNode(&left, data, n_l, first, level+1);
    SplitNode(&right, data, n_r, first_r, level+1);

    return;
}

struct Node *BuildTree(double data[][D]) {
    
    struct Node root;

    return &root;
}

double TestPoint(struct Node *root, double *data) {
    
    struct Node *node = root;
    while (node->left != NULL && node->right != NULL) {
        ind = node->index;
        if (data[ind] < node->threshold)
            node = node->left;
        else
            node = node->right;
    }

    return node->label;
}

void TreeFree(struct Node *root) {
    if (root->right)

    return;
}

/*
struct node *CountSplit(float data[][d], int n, int a, node *) {

 * data = node data
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
 *
    
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

*/




int main(int argc, char *argv[]) { 

    double data[20][D] = {{3, 2, 1},
        {4, 1, 0},
        {3, 1, 0},
        {5, 3, 1},
        {1, 6, 1},
        {8, 2, 0},
        {2, 4, 1},
        {5, 2, 0},
        {9, 1, 0},
        {5, 4, 0},
        {2, 1, 0},
        {1, 1, 1},
        {3, 1, 1},
        {2, 1, 0},
        {8, 3, 0},
        {5, 1, 0},
        {6, 4, 0},
        {9, 1, 1},
        {6, 0, 1},
        {11, 11, 0}};

    Sort(data, 0, 19, 1);
    int i, j;
    for (i = 0; i < 20; ++i) {
        for (j = 0; j < D; ++j)
            printf("%f ", data[i][j]);
        printf("\n");
    }

    

    return 0;
}






