#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

// Class kDTree
//private:
void kDTree::mergeSort(vector<vector<int>> &list, int d) {
    if (list.size() < 2) return;
    int mid = list.size()/2;
    vector<vector<int>> leftList(list.begin(), list.begin() + mid), 
                        rightList(list.begin() + mid, list.end());
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

    pR->left = buildTree(leftList, (d + 1) % k);
    pR->right = buildTree(rightList, (d + 1) % k);
    return pR;
}

//public:
string kDTree::inorderTraversal(kDTreeNode* root) const {
    if (!root) return "";
    stringstream ss;
    ss << inorderTraversal(root->left);
    ss << *root << " ";
    ss << inorderTraversal(root->right);
    return ss.str();
};
void kDTree::inorderTraversal() const{
    string res = this->inorderTraversal(this->root);
    res.pop_back();
    cout << res;
};

string kDTree::preorderTraversal(kDTreeNode* root) const {
    if (!root) return "";
    stringstream ss;
    ss << *root << " ";
    ss << preorderTraversal(root->left);
    ss << preorderTraversal(root->right);
    return ss.str();
};
void kDTree::preorderTraversal() const{
    string res = this->preorderTraversal(this->root);
    res.pop_back();
    cout << res;
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
    string res = this->postorderTraversal(this->root);
    res.pop_back();
    cout << res;
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
    this->k = pointList[0].size();
    this->root = buildTree(list, 0);
    return;
};

// count square of distance of 2 vector
long long int sqr_distance(const vector<int> &v1, const vector<int> &v2) {
    long long int res = 0;
    int len = v1.size();
    for (int i = 0; i < len; i++){
        int tmp = v1[i] - v2[i];
        res += pow(tmp, 2);
    }
    return res;
};

// to find closest from n1(traver 1st) and n2(traverse after n1) to target
kDTreeNode* closest_toTarget(kDTreeNode* n1, kDTreeNode* n2, const vector<int> &target){
    if (!n1) return n2;
    if (!n2) return n1;

    long long int   d1 = sqr_distance(n1->data, target), 
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
    }
    mergesortList(list,target);
    if (k < (int)list.size()) list.resize(k);
    return pR;
};

void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList){
    knearestNeighbour(this->root, target, k, 0, bestList);
    return;
};

kDTreeNode* kDTree::buildTree_v2(vector<vector<int>> &v_X, int d, const int k_dim) {
    if (v_X.empty()) return nullptr;
    int mid = ((int)v_X.size() - 1)/2;
    mergeSort(v_X, d);
    int lb = v_X[mid].back();
    v_X[mid].pop_back();
    kDTreeNode* pR = new kDTreeNode(v_X[mid], lb);

    // split pointList to left and right

    vector<vector<int>> leftList, rightList;
    for (int z = 0; z < mid; z++) {
        leftList.push_back(v_X[z]);
    }
    for(size_t z = mid + 1; z < v_X.size(); z++) {
        rightList.push_back(v_X[z]);
    }

    // recursion
    pR->left = buildTree_v2(leftList, (d + 1) % k_dim, k_dim);
    pR->right = buildTree_v2(rightList, (d + 1) % k_dim, k_dim);
    return pR;
}

void kDTree::buildTree_v2(const vector<vector<int>> &v_X, const vector<vector<int>> &v_y) {
    if (this->root) clear(this->root);
    this->count = v_X.size();
    if (count == 0) {
        this->root = nullptr;
        return;
    }
    vector<vector<int>> list_X(v_X.begin(),v_X.end());

    for (size_t i = 0; i < v_X.size(); i++) {
        list_X[i].push_back(v_y[i].front());
    }
    int k_dim = list_X.front().size() - 1;
    this->k = k_dim;
    this->root = buildTree_v2(list_X, 0, k_dim);
    return;
}
// end Class kDTree

// Class kNN

void kNN::fit(Dataset &X_train, Dataset &y_train) {
    vector<vector<int>> v_y, v_X; //to build X_tree and y_tree
    
    // transform Dataset y_train -> vector<vector<int>>, pre-condition of buildtree function
    while (y_train.data.size()){
        vector<int> tmp { make_move_iterator(begin(y_train.data.front())),
                          make_move_iterator(end(y_train.data.front())) };
        y_train.data.pop_front();
        v_y.push_back(tmp);
    }

    // transform Dataset X_train -> vector<vector<int>>, pre-condition of buildtree function
    while (X_train.data.size()) {
        // vector<int> tmp(X_train.data.front().begin(),X_train.data.front().end());
        vector<int> tmp { make_move_iterator(begin(X_train.data.front())),
                          make_move_iterator(end(X_train.data.front())) };
        X_train.data.pop_front();
        v_X.push_back(tmp);
    }

    // build tree
    this->X_tree.buildTree_v2(v_X,v_y);
}

Dataset kNN::predict(Dataset &X_test) {
    int len = X_test.data.size();
    Dataset y_pred;
    y_pred.columnName.push_back("label");
    
    //run row by row X_test
    for (int i = 0; i < len; i++) {
        vector<kDTreeNode*> bestList; // vector knearest neighbour

        vector<int> target { make_move_iterator(begin(X_test.data.front())),
                             make_move_iterator(end(X_test.data.front()))   };
        X_tree.kNearestNeighbour(target, this->k, bestList);
        // vector save label of bestList
        vector<int> findMode_label;
        for (int j = 0; j < this->k; j++) {
            findMode_label.push_back(bestList[j]->label);
        }
        // reference array to count freq of bestList label
        int ref_arr[10] = { 0 };

        for (int j = 0; j < this->k; j++) {
            ref_arr[findMode_label[j]]++;
        }
        // find Mode label of bestList
        int max = 0, idx = 0;
        for (int j = 0; j < 10; j++) {
            if (ref_arr[j] > max) {
                idx = j;
                max = ref_arr[j];
            }
        }
        // push predict value to y_pred
        list<int> temp;
        temp.push_back(idx);
        y_pred.data.push_back(temp);
        // pop the first element out of list
        X_test.data.pop_front();
    }
    return y_pred;
}

double kNN::score(const Dataset &y_test, const Dataset &y_pred) {
    int S = y_test.data.size();
    Dataset y_T = y_test, y_P = y_pred;
    int n_CorrectImg = 0;
    for (int i = 0; i < S; i++) {
        if (y_T.data.front().front() == y_P.data.front().front()) {
            n_CorrectImg++;
        }
        y_T.data.pop_front();
        y_P.data.pop_front();
    }
    double accuracy = (double)n_CorrectImg / (double)S;
    return accuracy;
}
// end Class Knn