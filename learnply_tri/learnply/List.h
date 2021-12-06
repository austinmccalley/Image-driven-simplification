#ifndef __LIST__H__
#define __LIST__H__

struct Node {
	void* value;
	struct Node* next;
};

struct LinkedList {
	Node* head;
	Node* tail;
	int size;
};

void createList(LinkedList *list);
void addToList(LinkedList *l, void *v);
void removeFromList(LinkedList *l, void *v);
int containsVal(LinkedList *l, void *v);
void destroyList(LinkedList *l);

#endif  //!__LIST__H__