#include <stdlib.h>
#include <stdio.h>
#include "tree.h"


void MergeSort(double **data, int first, int last, int a) {

    if (first >= last)
        return;

    int mid = (last-first)/2;
    MergeSort(data, first, first+mid, a);
    MergeSort(data, first+mid+1, last, a);

    Merge(data, first, first+mid, last, a);
    return;
}


void Merge(double **data, int first, int mid, int last, int a) {
    //TODO

    return;
}


void Sort(double **data, int first, int last, int a) {

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


void Partition(double **data, int first, int last, int a, int ends[]) {

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

