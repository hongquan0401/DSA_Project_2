#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

// Class kDTree
void kDTree::inorderTraversal(kDTreeNode* root, int count) const {
    if (!root) return;
    inorderTraversal(root->left, --count);
    if (count <= 0) cout << *root;
    else cout << *root << " ";
    inorderTraversal(root->right, --count);
};
void kDTree::inorderTraversal() const{
    this->inorderTraversal(this->root,this->count);
};

void kDTree::preorderTraversal(kDTreeNode* root, int count) const {
    if (!root) return;
    if (count <= 0) cout << *root;
    else cout << *root << " ";
    preorderTraversal(root->left, --count);
    preorderTraversal(root->right, --count);
};
void kDTree::preorderTraversal() const{
    this->preorderTraversal(this->root,this->count);
};

void kDTree::postorderTraversal(kDTreeNode* root, int count) const {
    if (!root) return;
    postorderTraversal(root->left, --count);
    postorderTraversal(root->right, --count);
    if (count <= 0) cout << *root;
    else cout << *root << " ";
};
void kDTree::postorderTraversal() const{
    this->postorderTraversal(this->root,this->count);
};

int kDTree::height(kDTreeNode* root) const{
    if (!root) return 0;
    else return max(height(root->left),height(root->right)) + 1;
}

int kDTree::leafCount(kDTreeNode* root) const{
    if(!root) return 0;
    if (!root->left && !root->right) return 1;
    return leafCount(root->left) + leafCount(root->right);
}
// end Class kDTree