#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

// Class kDTree
//private:
void kDTree::mergeSort(vector<vector<int>> &list, int d) {
    if (list.size() < 2) return;
    int mid = list.size()/2;
    vector<vector<int>> leftList(list.begin(), list.begin()+mid), 
                        rightList(list.begin()+mid, list.end());
    mergeSort(leftList, d);
    mergeSort(rightList, d);
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
    mergeSort(pointList, d);
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

    pR->left = buildTree(leftList, (d+1)%k);
    pR->right = buildTree(rightList, (d+1)%k);
    return pR;
}

//public:
string kDTree::inorderTraversal(kDTreeNode* root
) const {
    if (!root) return "";
    stringstream ss;
    ss << inorderTraversal(root->left);
    ss << *root << " ";
    ss << inorderTraversal(root->right);
    return ss.str();
};
void kDTree::inorderTraversal() const{
    cout << this->inorderTraversal(this->root);
};

string kDTree::preorderTraversal(kDTreeNode* root) const {
    if (!root) return "";
    stringstream ss;
    ss << *root << " ";
    ss << preorderTraversal(root->left);
    ss << preorderTraversal(root->right);
    string res = ss.str();
    res.resize(res.size() - 1);
    return res;
};
void kDTree::preorderTraversal() const{
    cout << this->preorderTraversal(this->root) << "o";
};

string kDTree::postorderTraversal(kDTreeNode* root) const {
    if (!root) return "";
    stringstream ss;
    ss << postorderTraversal(root->left);
    ss << postorderTraversal(root->right);
    ss << *root << " ";
    return ss.str();
};
void kDTree::postorderTraversal() const{
    cout << this->postorderTraversal(this->root);
};

int kDTree::height(kDTreeNode* root) const{
    if (!root) return 0;
    else return max(height(root->left), height(root->right)) + 1;
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
        return;
    }
    // kDTreeNode* tmp = removePoint(point);
    // if (!tmp && this->count > 1) this->count--;
    // else if (!tmp && this->count <= 1) this->count = 0;
    // else this->count--;
    this->root = removePoint(this->root,point,0);
    this->count--;
    return;
};

bool kDTree::search(const vector<int> &point){
    return search(this->root, point, 0);
};

void kDTree::buildTree(const vector<vector<int>> &pointList){
    this->count = pointList.size();
    
    if (this->root) clear(this->root);
    if (count == 0) {
        this->root = nullptr;
        return;
    }
    
    vector<vector<int>> list(pointList.begin(), pointList.end());
    this->root = buildTree(list, 0);
    return;
};

// count square of distance of 2 vector
long int sqr_distance(const vector<int> &v1, const vector<int> &v2) {
    long int res = 0;
    for (size_t i = 0; i < v1.size(); i++){
        res += pow((v1[i] - v2[i]), 2); 
    }
    return res;
};

// to find closest from n1(traver 1st) and n2(traverse after n1) to target
kDTreeNode* closest_toTarget(kDTreeNode* n1, kDTreeNode* n2, const vector<int> &target){
    if (!n1) return n2;
    if (!n2) return n1;

    long int    d1 = sqr_distance(n1->data, target), 
                d2 = sqr_distance(n2->data, target);
    if (d1 < d2) return n1;
    else return n2;
}

kDTreeNode* kDTree::nearestNeighbour(kDTreeNode* pR, const vector<int> &target, int d){
    if (!pR) return nullptr;
    kDTreeNode* nextBranch = nullptr, * otherBranch = nullptr;
    // to determine next branch to traverse
    if (target[d] < pR->data[d]) {
        nextBranch = pR->left;
        otherBranch = pR->right;
    }
    else {
        nextBranch = pR->right;
        otherBranch = pR->left;
    }
    // recursion
    kDTreeNode* temp = nearestNeighbour(nextBranch, target, (d + 1) % this->k);
    kDTreeNode* best = closest_toTarget(pR, temp, target);
    // check if R>= r then change to otherBranch
    long int R_sqr = sqr_distance(target, best->data), 
             r_sqr = pow(target[d] - pR->data[d], 2);
    if (R_sqr >= r_sqr) {
        temp = nearestNeighbour(otherBranch, target, (d + 1) % this->k);
        best = closest_toTarget(best, temp, target); // update best 
    }
    return best;
};

void kDTree::nearestNeighbour(const vector<int> &target, kDTreeNode *&best){
    best = nearestNeighbour(this->root, target, 0);
    return;
};

void mergesortList(vector<kDTreeNode*> &list, const vector<int> &target) {
    if (list.size() < 2) return;
    int mid = list.size() / 2;
    vector<kDTreeNode*> leftList(list.begin(), list.begin()+mid), 
                        rightList(list.begin()+mid, list.end());
    mergesortList(leftList, target);
    mergesortList(rightList, target);
    int nL = leftList.size(),
        nR = rightList.size();
    int i = 0, j = 0, k = 0; //    i for list, j for leftList, k for rightList
    while (j < nL && k < nR) {
        if (sqr_distance(leftList[j]->data,target) < sqr_distance(rightList[k]->data,target))
            list[i++] = leftList[j++];
        else list[i++] = rightList[k++];
    }
    while (j < nL) list[i++] = leftList[j++];
    while (k <nR) list[i++] = rightList[k++];
};

kDTreeNode* kDTree::knearestNeighbour(kDTreeNode* pR, const vector<int> &target, int k, int d,
                                      vector<kDTreeNode*> &list)
{
    if (!pR) return nullptr;
    kDTreeNode* nextBranch = nullptr, * otherBranch = nullptr;
    // to determine next branch to traverse
    if (target[d] < pR->data[d]) {
        nextBranch = pR->left;
        otherBranch = pR->right;
    }
    else {
        nextBranch = pR->right;
        otherBranch = pR->left;
    }
    // recursion
    kDTreeNode* temp = knearestNeighbour(nextBranch, target, k, (d + 1) % this->k, list);
    // push travered node
    list.push_back(pR);
    kDTreeNode* best = closest_toTarget(pR, temp, target);
    
    // check if R>= r then change to otherBranch
    long int R_sqr = sqr_distance(target, best->data), 
             r_sqr = pow(target[d] - pR->data[d], 2);
    if (R_sqr >= r_sqr) {
        temp = knearestNeighbour(otherBranch, target, k, (d + 1) % this->k, list);
        //best = closest_toTarget(best, temp, target); // update best 
    }
    mergesortList(list,target);
    if (k < (int)list.size()) list.resize(k);
    return pR;
};

void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList){
    knearestNeighbour(this->root, target, k, 0, bestList);
    return;
};
// end Class kDTree