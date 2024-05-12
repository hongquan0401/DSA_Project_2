#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

// Class kDTree
//private:
void kDTree::mergeSort(vector<vector<int>> &list, int d) {
    if (list.size() < 2) return;
    int mid = list.size()/2;
    vector<vector<int>> leftList(list.begin(),list.begin()+mid),
                        rightList(list.begin()+mid,list.end());
    mergeSort(leftList,d);
    mergeSort(rightList,d);
    int nL = leftList.size(),
        nR = rightList.size();
    int i = 0, j = 0, k = 0; // i for list, j for leftList, k for rightList
    while (j < nL && k < nR) { // insert while both leftList and rightList not empty
        if (leftList[j][d] < rightList[k][d]) list[i++] = leftList[j++];
        else list[i++] = rightList[k++];
    }
    while (j < nL) list[i++] = leftList[j++];
    while (k < nR) list[i++] = rightList[k++];
};
kDTreeNode* kDTree::buildTree(vector<vector<int>> &pointList, int d) {
    if (pointList.empty()) return nullptr;
    int mid = (pointList.size() - 1)/2;
    mergeSort(pointList,d);
    kDTreeNode* pR = new kDTreeNode(pointList[mid]);

    // split pointList to left and right

    vector<vector<int>> leftList, rightList;
    for (int z = 0; z < mid; z++) {
        leftList.push_back(pointList[z]);
    }
    for(size_t z = mid + 1; z < pointList.size(); z++) {
        rightList.push_back(pointList[z]);
    }

    // recursion

    pR->left = buildTree(leftList,(d+1)%k);
    pR->right = buildTree(rightList,(d+1)%k);
    return pR;
}

//public:
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
    if (!search(point)) {
        cout << "nothing";
        return;
    }
    kDTreeNode* tmp = removePoint(point);
    if (!tmp && this->count > 1) this->count--;
    else if (!tmp && this->count <= 1) this->count = 0;
    else this->count--;
    cout << this->count << " ";
    return;
};

bool kDTree::search(const vector<int> &point){
    return search(this->root,point,0);
};

void kDTree::buildTree(const vector<vector<int>> &pointList){
    this->count = pointList.size();
    
    if (this->root) clear(this->root);
    if (count == 0) {
        this->root = nullptr;
        return;
    }
    
    vector<vector<int>> list(pointList.begin(),pointList.end());
    this->root = buildTree(list,0);
    return;
};
void kDTree::nearestNeighbour(const vector<int> &target, kDTreeNode *&best){
    return;
};
void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList){
    return;
};
// end Class kDTree