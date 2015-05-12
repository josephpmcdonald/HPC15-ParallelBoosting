#define GINI(pos, n) (2.*(pos)/(n)*((n)-(pos))/(n))
#define ENTROPY(pos, n) ((pos)/(n)*log((pos)/(n)) + ((n)-(pos))/(n)*log(((n)-(pos))/(n)))
#define WGINI(pos, tot) (2.*(pos)/(tot)*((tot)-(pos))/(tot))
#define D 3


struct TreeNode {
    struct TreeNode *parent;
    struct TreeNode *right;
    struct TreeNode *left;
    int index;
    double threshold;
    double label;
};

typedef struct TreeNode Node;

//sort.c
void MergeSort(double **data, int first, int last, int a);
void Merge(double **data, int first, int mid, int last, int a);
void Sort(double **data, int first, int last, int a);
void Partition(double **data, int first, int last, int a, int ends[]);

//tree.c
int BestSplit(double **data, int n, int first, int col, int pos, double *impurity);
void SplitNode(Node *node, double **data, int n, int first, int level);
void BuildTree(Node *root, double **data, int n);
double TestPoint(Node *root, double *data);
void TreeFree(Node *node);
void TreePrint(Node *node, int level);


//data.c
double **MNIST17();




