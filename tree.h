#define GINI(pos, n) (2.*(pos)/(n)*((n)-(pos))/(n))
#define ENTROPY(pos, n) ((pos)/(n)*log((pos)/(n)) + ((n)-(pos))/(n)*log(((n)-(pos))/(n)))
#define D 3


struct Node {
    struct Node *parent;
    struct Node *right;
    struct Node *left;
    int index;
    double threshold;
    double label;
};

//tree.c
void Sort(double data[][D], int first, int last, int a);
void Partition(double data[][D], int first, int last, int a, int ends[]);
int BestSplit(double data[][D], int n, int first, int col, int pos, double *impurity);
void SplitNode(struct Node *node, double data[][D], int n, int first, int level);

//data.c
double **MNIST17();




