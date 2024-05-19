#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    int label;
    kDTreeNode(vector<int> data, int label = 0, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
        this->label = label;
    }

    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for (size_t i = 0; i < node.data.size(); i++)
        {
            os << node.data[i];
            if (i != node.data.size() - 1)
            {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;
    int count;

    void clear(kDTreeNode* root) {
        if (!root) return;
        else {
            if (root->left) clear(root->left);
            if (root->right) clear(root->right);
            delete root;
            return;
        }
    };

    kDTreeNode* deepCopyNode(const kDTreeNode* pR) {
        if (!pR) return nullptr;
        kDTreeNode* pNew = new kDTreeNode(pR->data);
        pNew->left = deepCopyNode(pR->left);
        pNew->right = deepCopyNode(pR->right);
        return pNew;
    };
    string inorderTraversal(kDTreeNode* root) const;
    string preorderTraversal(kDTreeNode* root) const;
    string postorderTraversal(kDTreeNode* root) const;
    int height(kDTreeNode* root) const;
    int leafCount(kDTreeNode* root) const;
    void runDim(int& dim) { if (++dim >= k) dim = 0; }
    bool samePoint(const vector<int> &p1, const vector<int> &p2) {
        for (int i = 0; i < k; i++) {
            if (p1[i] != p2[i]) return false;
        }
        return true;
    }
    bool search(kDTreeNode* pR, const vector<int> &point, int dim) {
        if (!pR) return false;
        if (samePoint(pR->data, point)) return true;
        dim %= this->k;
        if (point[dim] < pR->data[dim]) { 
            return search(pR->left, point, ++dim);
        }
        else {
            return search(pR->right, point, ++dim);
        }
    }
    // find min point in particular dim value
    kDTreeNode* minPoint(kDTreeNode* x, kDTreeNode* y, kDTreeNode* z, int dim){
        kDTreeNode* res = x;
        if (y && y->data[dim] < res->data[dim]) res = y;
        if (z && z->data[dim] < res->data[dim]) res = z;
        return res;
    }
    // find min Node to replace remove Node
    kDTreeNode* findMin(kDTreeNode* root, int alpha, int d){
        if (!root) return nullptr;
        int cd = d % this->k;
        if (alpha == cd) {
            if (!root->left) return root;
            return findMin(root->left, alpha, d + 1);
        }
        return minPoint(root, findMin(root->left, alpha, d + 1), findMin(root->right, alpha, d + 1), alpha);
    }
    kDTreeNode* removePoint(kDTreeNode* pR, const vector<int> &point, int d) {
        if (!pR) return nullptr;
        int cd = d % this->k;
        // if root->data = point
        if (samePoint(pR->data, point)) {
            if (pR->right) {
                kDTreeNode* min = findMin(pR->right, cd, cd + 1);
                // copyPoint(pR->data, min->data);
                pR->data = min->data;
                pR->right = removePoint(pR->right, min->data, d + 1);
            }
            else if (pR->left){
                kDTreeNode* min = findMin(pR->left, cd, cd + 1);
                // copyPoint(pR->data, min->data);
                pR->data = min->data;
                pR->right = removePoint(pR->left, min->data, d + 1);
                pR->left = nullptr;
            }
            else {
                delete pR;
                return nullptr;
            }
            return pR;
        }
        if (point[cd] < pR->data[cd])
            pR->left = removePoint(pR->left, point, d + 1);
        else 
            pR->right = removePoint(pR->right, point, d + 1);
        return pR;
    };

    void mergeSort(vector<vector<int>> &list, int d);
    kDTreeNode* buildTree(vector<vector<int>> &pointList, int d);
    kDTreeNode* buildTree_v2(vector<vector<int>> &v_X, int d, const int k_dim);

    kDTreeNode* nearestNeighbour(kDTreeNode* pR, const vector<int> &target, int d);
    kDTreeNode* knearestNeighbour(kDTreeNode* pR, const vector<int> &target, int k, int d, 
                                        vector<kDTreeNode*> &list);
public:
    kDTree(int k = 2): k(k), root(nullptr), count(0) {};
    ~kDTree() {
        this->clear(this->root);
    };

    const kDTree &operator=(const kDTree &other) {
        this->count = other.count;
        this->k = other.k;
        this->root = deepCopyNode(other.root);
        return *this;
    };

    kDTree(const kDTree &other){
        this->count = other.count;
        this->k = other.k;
        this->root = deepCopyNode(other.root);
    };

    void inorderTraversal() const;
    void preorderTraversal() const;
    void postorderTraversal() const;
    int height() const { return height(this->root); };
    int nodeCount() const { return count; };
    int leafCount() const { return leafCount(this->root); };

    void insert(const vector<int> &point);
    void remove(const vector<int> &point);
    bool search(const vector<int> &point);
    void buildTree(const vector<vector<int>> &pointList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);

    kDTreeNode* findMin(kDTreeNode* pR, const int dim) { return findMin(pR, dim, 0); };
    kDTreeNode* findMin(const int dim) { return findMin(this->root, dim, 0); };
    void buildTree_v2(const vector<vector<int>> &v_X, const vector<vector<int>> &v_y);
};

class kNN
{
private:
    int k;
    kDTree X_tree;
public:
    kNN(int k = 5): k(k), X_tree(kDTree()) {}
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed
