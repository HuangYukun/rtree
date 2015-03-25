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

}

//when calling this function we should pass in the root node
RTNode * RTree::choose_leaf(const vector<int>& coordinate, RTNode* N)//recursive call
{
	// RTNode *N = this->root;
	if (N->entry_num<2||(N->entries[0]).get_ptr()==NULL)//if N is a leaf
	{
		return N;
	}
	else
	{
		vector<int> enlargement;
		for (int i = 0; i < N->entry_num; i++){
			BoundingBox bb = BoundingBox((N->entries[i]).get_mbr());
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
		//print min_index to debug
		choose_leaf(coordinate, (N->entries[min_index.back()]).get_ptr());
	}

}

void RTree::adjust_tree(RTNode* L, RTNode* LL)
{
	RTNode* N = L;
	RTNode* NN = LL;
	//check N is the root
	if (N->parent==NULL){
		return;
	}
	else
	{
		RTNode* P = N->parent;
		for (int i = 0; i < P->entry_num; i++){
			if(P->entries[i].get_ptr()==N){
				for (int j = 0; j < N->entry_num; j++){
					BoundingBox mbr(P->entries[i].get_mbr());
					mbr.group_with((N->entries[j]).get_mbr());
					P->entries[i].set_mbr(mbr);
				}
			}
		}
		if (NN!=NULL){
			Entry enn = Entry();
			enn.set_ptr(NN);
			BoundingBox bnn = get_mbr(NN->entries, NN->entry_num);
			enn.set_mbr(bnn);
			if (P->entry_num < P->size){
				//add enn to P
				P->entries[P->entry_num] = enn;
				P->entry_num++;
			}
			else{
				//AT4
				//firstly overflow P, then split it
				P->entries[P->entry_num] = enn;
				P->entry_num++;
				vector<RTNode*> Ps = split_node(P);
				N = Ps.front();
				NN = Ps.back();
			}
		}
	}
}

//linear-cost algorithm
vector<RTNode*> RTree::split_node(RTNode* node){
	vector<int> separation;
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
		separation.push_back(abs(*highest_low_side - *lowest_high_side)/(*highest_high_side - *lowest_low_side));
	}
	vector<int>::iterator biggest_normalized_separation = max_element(separation.begin(),separation.end());
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
	RTNode * LL = new RTNode(node->level, node->size);
	LL->entries[0] = node->entries[lhs_entry_index];
	LL->entry_num++;
	for (int i = lhs_entry_index; i < node->size; i++){
		node->entries[i] = node->entries[i+1];
	}
	node->entry_num--;
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
