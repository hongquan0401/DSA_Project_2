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
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }

    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for (int i = 0; i < node.data.size(); i++)
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
        if (root->left) clear(root->left);
        if (root->right) clear(root->right);
        delete root;
    }

    kDTreeNode* deepCopyNode(const kDTreeNode* pR) {
        if (!pR) return nullptr;
        kDTreeNode* pNew = new kDTreeNode(pR->data);
        pNew->left = deepCopyNode(pR->left);
        pNew->right = deepCopyNode(pR->right);
        return pNew;
    };
    void inorderTraversal(kDTreeNode* root, int &count) const;
    void preorderTraversal(kDTreeNode* root, int &count) const;
    void postorderTraversal(kDTreeNode* root, int &count) const;
    int height(kDTreeNode* root) const;
    int leafCount(kDTreeNode* root) const;
public:
    kDTree(int k = 2): k(k), root(nullptr), count(0) {};
    ~kDTree() {
        if (this->root) {
            clear(this->root);
            count = 0;
        }
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
};

class kNN
{
private:
    int k;

public:
    kNN(int k = 5);
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed
