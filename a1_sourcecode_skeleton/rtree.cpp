/* Implementations of R tree */
#include <cmath>
#include <algorithm>
#include "rtree.h"


const double EPSILON = 1E-10;

RTree::RTree(int entry_num)
{
	max_entry_num = entry_num;
	dimension = 2;//by default
	root = new RTNode(0, entry_num);
}

RTree::RTree(int entry_num, int dim)
{
	max_entry_num = entry_num;
	dimension = dim;//by default
	root = new RTNode(0, entry_num);
}

RTree::~RTree()
{
	delete root;
	root = NULL;
}


bool RTree::insert(const vector<int>& coordinate, int rid)
{
	if (coordinate.size() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}

	/***
	ADD YOUR CODE HERE
	****/

	//need to handle duplicated case
	Entry qpEmptyEntry;
	if (query_point(coordinate, qpEmptyEntry)==true){
		return false;
	}
	//can be done as first making a query
	//start of insert
	RTNode * L = new RTNode(0, coordinate.size());
	L = choose_leaf(coordinate, root);
	RTNode * LL = NULL;
	BoundingBox bb = BoundingBox(coordinate, coordinate);
	Entry entry_to_insert = Entry(bb, rid);
	if (L->entry_num < L->size){
		L->entries[L->entry_num] = Entry(bb, rid);
		L->entry_num++;
	}
	else{
		vector<RTNode*> splitted_nodes = split_node(L, entry_to_insert);
		L = splitted_nodes.front();
		LL = splitted_nodes.back();
	}
	vector<RTNode*> roots = adjust_tree(L, LL);

	if (roots.back()!=NULL){
		//create a new root with two entries pointing to L&LL
		RTNode* new_root = new RTNode(roots.back()->level+1, roots.back()->size);
		roots.front()->parent = new_root;
		roots.back()->parent = new_root;
		new_root->entries[0] = Entry();
		new_root->entries[0].set_mbr(get_mbr(roots.front()->entries, roots.front()->entry_num));
		new_root->entries[0].set_ptr(roots.front());
		new_root->entries[1] = Entry();
		new_root->entries[1].set_mbr(get_mbr(roots.back()->entries, roots.back()->entry_num));
		new_root->entries[1].set_ptr(roots.back());
		root = new_root;
		root->entry_num = 2;
	}

	return true;

}

//when calling this function we should pass in the root node
RTNode * RTree::choose_leaf(const vector<int>& coordinate, RTNode* N)//recursive call
{
	// RTNode *N = this->root;
	if (N->level==0)//if N is a leaf
	{
		return N;
	}
	else
	{
		vector<int> enlargement;
		for (int i = 0; i < N->entry_num; i++){
			BoundingBox bb = BoundingBox((N->entries[i]).get_mbr());
			//basically the entries are points with not range
			BoundingBox box_new = BoundingBox(coordinate,	coordinate);
			bb.group_with(box_new);
			enlargement.push_back(bb.get_area() - (N->entries[i]).get_mbr().get_area());
		}
		//tie-breaking
		int min = enlargement.front();
		vector<int> min_index;
		for (vector<int>::iterator it = enlargement.begin() ; it != enlargement.end(); ++it){
			if (*it < min){
				min_index.clear();
				min_index.push_back(it - enlargement.begin());
			}
			else if(*it == min){
				min_index.push_back(it - enlargement.begin());
			}
		}
		for (vector<int>::iterator it = min_index.begin()+1 ; it != min_index.end(); ++it){
			BoundingBox b1 = (N->entries[*(it-1)]).get_mbr();
			BoundingBox b2 = (N->entries[*it]).get_mbr();
			if (tie_breaking(b1, b2) == true){
				min_index[*it] = min_index[*(it-1)];
			}
		}
		//the return cost me an hour
		return choose_leaf(coordinate, (N->entries[min_index.back()]).get_ptr());
	}

}

//return the single root or splitted roots
vector<RTNode*> RTree::adjust_tree(RTNode* L, RTNode* LL) //LL can be NULL
{
	//assume the count starts from 1
	// RTNode* N = new RTNode(L->level, L->size);
	RTNode* N = L;
	RTNode* NN = LL; //check NULL
	//check N is the root
	//parent pointer need taken care
	if (N->parent==NULL){
		vector<RTNode*> roots;
		roots.push_back(N);
		//some problems with this
		roots.push_back(NN);
			// cout << "N is " << N->level <<endl;
			// cout << "NN is " << NN->level <<endl;
		return roots;
	}
	else
	{
		RTNode* P = new RTNode(N->level+1, N->size);
		P = N->parent;
		for (int i = 0; i < P->entry_num; i++){
			if(P->entries[i].get_ptr()==N){
				for (int j = 0; j < N->entry_num; j++){
					BoundingBox mbr(P->entries[i].get_mbr());
					mbr.group_with((N->entries[j]).get_mbr());
					P->entries[i].set_mbr(mbr);
				}
			}
		}
		// AT4 + AT5, not a MUST
		if (NN!=NULL){
			Entry enn = Entry();
			enn.set_ptr(NN);
			BoundingBox bnn = get_mbr(NN->entries, NN->entry_num);
			enn.set_mbr(bnn);
			if (P->entry_num < P->size){
				//add enn to P with room
				P->entries[P->entry_num] = Entry();
				P->entries[P->entry_num] = enn;
				P->entry_num++;
				N = P;
				return adjust_tree(N, NULL);

			}
			else{
				//AT4
				vector<RTNode*> Ps = split_node(P, enn);
				N = Ps.front();
				NN = Ps.back();
					// cout << "N " << N->level <<endl;
					// cout << "NN " << NN->level <<endl;
				return adjust_tree(N, NN);
			}
		}
		//no NN is passed in
		else
		{
			N = P;
			return adjust_tree(N, NULL);
		}

	}
}

//linear-cost algorithm
//change it to take in new entry e
vector<RTNode*> RTree::split_node(RTNode* node, Entry entry_to_insert){
	vector<double> separation;
	for (int j = 0; j < node->entries[0].get_mbr().get_dim(); j++){
		vector<int> low_side;
		vector<int> high_side;
		for (int i = 0; i < node->entry_num; i++){
			low_side.push_back(node->entries[i].get_mbr().get_lowestValue_at(j));
			high_side.push_back(node->entries[i].get_mbr().get_highestValue_at(j));
		}
		vector<int>::iterator highest_low_side = max_element(low_side.begin(),low_side.end());
		vector<int>::iterator lowest_high_side = min_element(high_side.begin(),high_side.end());

		vector<int>::iterator lowest_low_side = min_element(low_side.begin(),low_side.end());
		vector<int>::iterator highest_high_side = max_element(high_side.begin(),high_side.end());
		// type casting, to double
		separation.push_back(abs(*highest_low_side - *lowest_high_side)/double((*highest_high_side - *lowest_low_side)));
	}
	vector<double>::iterator biggest_normalized_separation = max_element(separation.begin(),separation.end());
	//back trace to the two entries
	int dimension_index = biggest_normalized_separation - separation.begin();
	int highest_low_side = node->entries[0].get_mbr().get_lowestValue_at(dimension_index);
	int hls_entry_index = 0;
	int lowest_high_side = node->entries[0].get_mbr().get_highestValue_at(dimension_index);
	int lhs_entry_index = 0;
	for (int i = 0; i < node->entry_num; i++){
		int low_side = node->entries[i].get_mbr().get_lowestValue_at(dimension_index);
		if (low_side > highest_low_side){
			highest_low_side = low_side;
			hls_entry_index = i;
		}
		int high_side = node->entries[i].get_mbr().get_highestValue_at(dimension_index);
		if (high_side < lowest_high_side){
			lowest_high_side = high_side;
			lhs_entry_index = i;
		}
	}
	vector<RTNode*> Nodes;
	//modify constructor
	//new Node to carry splitted entries
	RTNode * LL = new RTNode(node->level, node->size);
	LL->entries[0] = Entry();
	LL->entries[0] = node->entries[lhs_entry_index];
	LL->entry_num++;
	//rearrange the entries sequence
	for (int i = lhs_entry_index; i < node->size-1; i++){
		node->entries[i] = node->entries[i+1];
	}
	node->entry_num--;

	//check where to input the new node
	//node
	int area_difference1;
	int area_difference2;
	BoundingBox bb1 = get_mbr(node->entries, node->entry_num);
	bb1.group_with(entry_to_insert.get_mbr());
	BoundingBox bb2 = get_mbr(LL->entries, LL->entry_num);
	bb2.group_with(entry_to_insert.get_mbr());
	area_difference1 = abs(get_mbr(node->entries, node->entry_num).get_area() - bb1.get_area());
	area_difference2 = abs(get_mbr(LL->entries, LL->entry_num).get_area() - bb2.get_area());
	if (area_difference1 < area_difference2){
		//insert into node
		node->entries[node->entry_num] = Entry();
		node->entries[node->entry_num] = entry_to_insert;
		node->entry_num++;
	}
	else{
		//insert into LL
		LL->entries[LL->entry_num] = Entry();
		LL->entries[LL->entry_num] = entry_to_insert;
		LL->entry_num++;
	}
	Nodes.push_back(node);
	Nodes.push_back(LL);
	return Nodes;
}

void RTree::query_range(const BoundingBox& mbr, int& result_count, int& node_travelled)
{
	if (mbr.get_dim() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}

	/***
	ADD YOUR CODE HERE
	****/
	RTNode * startNode = root;
	for (int i = 0; i < startNode->entry_num; i++){
		if (startNode->entries[i].get_mbr().is_intersected(mbr)==true){
			
		}
	}
}


bool RTree::query_point(const vector<int>& coordinate, Entry& result)
{
	if (coordinate.size() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}

	/***
	ADD YOUR CODE HERE
	****/

	RTNode * leaf = choose_leaf(coordinate, root);
	BoundingBox bb = BoundingBox(coordinate, coordinate);
	int index = -1;
	for (int i = 0; i < leaf->entry_num; i++){
		if (leaf->entries[i].get_mbr().is_equal(bb) == true){
			index = i;
		}
	}
	if (index >= 0 ){
		result = leaf->entries[index];
		return true;
	}
	return false;
}


/**********************************
 *
 * Please do not modify the codes below
 *
 **********************************/

//
// Calcuate the MBR of a set of entries, of size ``len''.
// Store the MBR in the first entry
//
BoundingBox RTree::get_mbr(Entry* entry_list, int len)
{
	BoundingBox mbr(entry_list[0].get_mbr());
	for (int i = 1; i < len; i++) {        
		mbr.group_with(entry_list[i].get_mbr());
	}
	return mbr;
}


/*********************************************************
  Return true means choose box1 for tie breaking.
  If the two boxes is the same, return true.
  This is to give a unified way of tie-breaking such that if your program is correct, then the result should be same, not influnced by any ties.
 *********************************************************/
bool RTree::tie_breaking(const BoundingBox& box1, const BoundingBox& box2)
{
	//for every dimension, try to break tie by the lowest value, then the highest
	for (int i = 0; i < box1.get_dim(); i++)
	{
		if (box1.get_lowestValue_at(i) != box2.get_lowestValue_at(i))
		{
			return box1.get_lowestValue_at(i) < box2.get_lowestValue_at(i);
		}
		else if (box1.get_highestValue_at(i) != box2.get_highestValue_at(i))
		{
			return box1.get_highestValue_at(i) > box2.get_highestValue_at(i);
		}
	}
	return true;
}


void RTree::stat(RTNode* node, int& record_cnt, int& node_cnt)
{
	if (node->level == 0) {
		record_cnt += node->entry_num;
		node_cnt++;
	}
	else {
		node_cnt++;
		for (int i = 0; i < node->entry_num; i++)
			stat((node->entries[i]).get_ptr(), record_cnt, node_cnt);
	}
}

void RTree::stat()
{
	int record_cnt = 0, node_cnt = 0;
	stat(root, record_cnt, node_cnt);
	cout << "Height of R-tree: " << root->level + 1 << endl;
	cout << "Number of nodes: " << node_cnt << endl;
	cout << "Number of records: " << record_cnt << endl;
	cout << "Dimension: " << dimension << endl;
}


void RTree::print_node(RTNode* node, int indent_level)
{
	BoundingBox mbr = get_mbr(node->entries, node->entry_num);

	char* indent = new char[4*indent_level+1];
	memset(indent, ' ', sizeof(char) * 4 * indent_level);
	indent[4*indent_level] = '\0';

	if (node->level == 0) {
		cout << indent << "Leaf node (level = " << node->level << ") mbr: (";
		for (int i = 0; i < mbr.get_dim(); i++)
		{
			cout << mbr.get_lowestValue_at(i) << " " << mbr.get_highestValue_at(i);
			if (i != mbr.get_dim() - 1)
			{
				cout << " ";
			}
		}
		cout << ")\n";
	}
	else {

		cout << indent << "Non leaf node (level = " << node->level << ") mbr: (";
		for (int i = 0; i < mbr.get_dim(); i++)
		{
			cout << mbr.get_lowestValue_at(i) << " " << mbr.get_highestValue_at(i);
			if (i != mbr.get_dim() - 1)
			{
				cout << " ";
			}
		}
		cout << ")\n";
	}

	Entry *copy = new Entry[node->entry_num];
	for (int i = 0; i < node->entry_num; i++) {
		copy[i] = node->entries[i];
	}

	for (int i = 0; i < node->entry_num; i++) {
		int index = 0; // pick next.
		for (int j = 1; j < node->entry_num - i; j++) {
			if (tie_breaking(copy[j].get_mbr(), copy[index].get_mbr())) {
				index = j;
			}
		}

		if (node->level == 0) {
			Entry& e = copy[index];
			cout << indent << "    Entry: <";
			for (int i = 0; i < e.get_mbr().get_dim(); i++)
			{
				cout << e.get_mbr().get_lowestValue_at(i) << ", ";
			}
			cout << e.get_rid() << ">\n";
		}
		else {
			print_node(copy[index].get_ptr(), indent_level+1);
		}
		// Move the output one to the rear.
		Entry tmp = copy[node->entry_num - i - 1];
		copy[node->entry_num - i - 1] = copy[index];
		copy[index] = tmp;

	}

	delete []indent;
	delete []copy;
}

void RTree::print_tree()
{
	if (root->entry_num == 0)
		cout << "The tree is empty now." << endl;
	else
		print_node(root, 0);
}
