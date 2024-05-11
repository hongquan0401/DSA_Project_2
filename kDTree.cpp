#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

// Class kDTree
void kDTree::inorderTraversal(kDTreeNode* root, int &count) const {
    if (!root) return;
    inorderTraversal(root->left, count);
    if (count <= 1) cout << *root;
    else {
            cout << *root << " ";
            count--;
    }
    inorderTraversal(root->right, count);
};
void kDTree::inorderTraversal() const{
    int tmp = this->count;
    this->inorderTraversal(this->root,tmp);
};

void kDTree::preorderTraversal(kDTreeNode* root, int &count) const {
    if (!root) return;
    if (count <= 1) cout << *root;
    else {
            cout << *root << " ";
            count--;
        }
    preorderTraversal(root->left, count);
    preorderTraversal(root->right, count);
};
void kDTree::preorderTraversal() const{
    int tmp = this->count;
    this->preorderTraversal(this->root,tmp);
};

void kDTree::postorderTraversal(kDTreeNode* root, int &count) const {
    if (!root) return;
    postorderTraversal(root->left, --count);
    postorderTraversal(root->right, --count);
    if (count <= 1) cout << *root;
    else {
            cout << *root << " ";
            count--;
    }
};
void kDTree::postorderTraversal() const{
    int tmp = this->count;
    this->postorderTraversal(this->root,tmp);
};

int kDTree::height(kDTreeNode* root) const{
    if (!root) return 0;
    else return max(height(root->left),height(root->right)) + 1;
}

int kDTree::leafCount(kDTreeNode* root) const{
    if(!root) return 0;
    if (!root->left && !root->right) return 1;
    return leafCount(root->left) + leafCount(root->right);
};

void kDTree::insert(const vector<int> &point) {
    kDTreeNode** pR = &this->root;
    int dim = 0;
    while (*pR) {
        pR = (point[dim] < (*pR)->data[dim]) ? &((*pR)->left) : &((*pR)->right);
        runDim(dim);
    } 
    *pR = new kDTreeNode(point);
    count++;
};

void kDTree::remove(const vector<int> &point){
    return;
};

bool kDTree::search(const vector<int> &point){
    return search(this->root,point,0);
};

void kDTree::buildTree(const vector<vector<int>> &pointList){
    return;
};
void kDTree::nearestNeighbour(const vector<int> &target, kDTreeNode *&best){
    return;
};
void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList){
    return;
};
// end Class kDTree