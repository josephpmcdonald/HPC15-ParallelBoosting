#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tree.h"


/* Notes: Use main to test functions.
 * 
 * TODO: Merge, TestPoint, tree pruning
 * 
 */


int BestSplit(double **data, int n, int first, int col, int pos, double *impurity) {

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

        while (i < n && (data[first+i][col] == threshold)) {
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


void SplitNode(Node *node, double **data, int n, int first, int level) {

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
    node->index = -1;

    //Get initial counts for positive/negative labels
    int i;
    int pos = 0;
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

    ///////////////TEST/////////////
    printf("pos = %d, neg = %d, lab = %f\n", pos, neg, node->label);


    ////////////////////////////////TEST//////////////
    //If branch is small or almost pure, make leaf
    printf("GINI: %f\n", GINI(pos, n));
    if (n < 6) {
        printf("small branch: %d points\n", n, level);
        return;
    }
    else if (level == 6) {
        printf("leaf node: level = 6\n");
        return;
    }
    else if (pos == 0 || neg == 0) {
        printf("pure node\n");
        return;
    }

    int col;
    int row; //best row to split at for particular column/feature
    int localrow; //first + localrow = row; receives BestSplit which returns integer in [-1, n-1]
    double threshold; //best threshold to split at for column/feature
    double impurity; //impurity for best split in feature/column
    int bestcol = -1; //feature with best split
    int bestrow = first+n-1; //best row to split for best feature
    double bestthresh; //threshold split for best feature (data[bestrow][bestcol])
    double Pmin = GINI(pos, n); //minimum impurity seen so far

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

//////////////////////////////////////////////////////////////////////
    printf("bestrow = %d, bestcol = %d\n", bestrow, bestcol);
    
    //If splitting doesn't improve purity (best split is at the end) stop
    if (bestrow == first+n-1) {
        printf("no improvement\n");
        return;
    }


    Sort(data, first, first+n-1, bestcol);

    //For feature, threshold with best impurity, save to node attributes
    node->index = bestcol;
    node->threshold = bestthresh;

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

    int first_r = bestrow+1;
    int n_l = first_r - first;
    int n_r = n - n_l;

    SplitNode(l, data, n_l, first, level+1);
    SplitNode(r, data, n_r, first_r, level+1);

    return;
}


void BuildTree(Node *root, double **data, int n) {

/* Wrapper function to initiate decision tree construction
 */

    SplitNode(root, data, n, 0, 0);
    return;
}


void TreeFree(Node *node) {

/* Recursively frees dynamically allocated trees
 */

    if (node == NULL)
        return;
    
    TreeFree(node->left);
    TreeFree(node->right);
    free(node);
    return;
}


double TestPoint(Node *root, double *data) {

    int ind;
    Node *node = root;
    while (node->left != NULL) {
        ind = node->index;
        if (data[ind] <= node->threshold)
            node = node->left;
        else
            node = node->right;
    }

    return node->label;
}


void TreePrint(Node *node, int level) {
    if (!node)
        return;

    TreePrint(node->left, level+1);
    printf("Node: level %d\n", level);
    printf("Index: %d\n", node->index);
    printf("Threshold: %f\n", node->threshold);
    printf("Class: %f", node->label);
    printf("\n");
    TreePrint(node->right, level+1);
    
    return;
}


int main(int argc, char *argv[]) { 

/*
    double **data17 = MNIST17();

    int i;
    for (i = 0; i < 10; ++i)
        printf("data17[%d]=%f\n", i, data17[i][28*28]);

    for (i = 0; i < 13007; ++i)
        free(data17[i]);
    free(data17);

    return 0;
*/


    double DATA[20][D] = {{3, 6, 1},
        {4, 1, -1},
        {3, 1, -1},
        {5, 7, 1},
        {1, 6, 1},
        {8, 2, 1},
        {2, 8, 1},
        {5, 2, -1},
        {9, 1, 1},
        {5, 4, -1},
        {2, 1, -1},
        {1, 6, 1},
        {3, 9, 1},
        {2, 1, -1},
        {8, 3, 1},
        {5, 1, -1},
        {6, 4, -1},
        {9, 5, 1},
        {6, 0, -1},
        {11, 11, 1}};

    double *data[20];
    int i;
    int j;
    for (i = 0; i < 20; ++i)
        data[i] = &DATA[i][0];

/*
    Sort(data, 0, 19, 1);

    double impurity;
    i = BestSplit(data, 20, 0, 1, 8, &impurity);
    printf("i = %d\n", i);
*/

    Node *root = malloc(sizeof(Node));
    root->right = NULL;
    root->left = NULL;
    root->parent = NULL;
    BuildTree(root, data, 20);

    for (i = 0; i < 20; ++i) {
        for (j = 0; j < D; ++j)
            printf("%f ", data[i][j]);
        printf("\n");
    }

    TreePrint(root, 0);

    double pt[] = {2, 4, -1};
    double *point = &pt[0];
    printf("Test: f(pt) = %f\n", TestPoint(root, point));

    TreeFree(root);

    return 0;
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



