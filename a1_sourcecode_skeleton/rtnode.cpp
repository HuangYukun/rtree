#include "rtnode.h"

//======================== Entry implementation =====================================================

Entry::Entry():mbr() {
	this->rid = -1;
	this->ptr = NULL;
}

Entry::Entry(const BoundingBox& thatMBR, const int rid):mbr(thatMBR) {
	this->rid = rid;
	this->ptr = NULL;
}

Entry::~Entry() {
	this->ptr = NULL;
}

const BoundingBox& Entry::get_mbr() const {
	return this->mbr;
}	


RTNode* Entry::get_ptr() const {
	return this->ptr;
}

int Entry::get_rid() const { 
	return this->rid;
}


void Entry::set_mbr(const BoundingBox& thatMBR) {
	this->mbr.set_boundingbox(thatMBR);//ATTENTION: should be deep copy...
}

void Entry::set_ptr(RTNode* ptr) {
	this->ptr = ptr;
}


void Entry::print() {
	this->mbr.print();
	cout << this->rid << endl;
	cout << this->ptr << endl;
}

//======================== RTNode implementation ==============================================

RTNode::RTNode(int lev, int s)
{
	entry_num = 0;
	entries = new Entry[s];
	//change it to s+1 for the sake of convenience
	// entries = new Entry[s+1];
	level = lev;
	size = s;
}

RTNode::RTNode(const RTNode& other)
{
	entries = new Entry[other.size];
	*this = other;
}

RTNode& RTNode::operator=(const RTNode& other)
{
	if (&other != this) {
		entry_num = other.entry_num;
		level = other.level;
		size = other.size;
		for (int i = 0; i < entry_num; i++)
			entries[i] = other.entries[i];
	}
	return *this;
}


RTNode::~RTNode()
{
	if (level != 0) {
		for (int i = 0; i < entry_num; i++) {
			delete entries[i].get_ptr();
			entries[i].set_ptr(NULL);
		}
	}
	delete []entries;
	entries = NULL;
}

// int RTNode::get_size()
// {
// 	return this->size;
// }
