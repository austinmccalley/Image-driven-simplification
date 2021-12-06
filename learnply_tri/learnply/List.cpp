#include "List.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;

void createList(LinkedList *l){
  l->head = NULL;
  l->tail = NULL;
  l->size = 0;
}

void addToList(LinkedList *l, void* v){
  Node *n = (Node*)malloc(sizeof(Node));
  n->value = v;
  n->next = NULL;
  if(l->size == 0){
    l->head = n;
    l->tail = n;
  }
  else{
    l->tail->next = n;
    l->tail = n;
  }
  l->size++;
}

void removeFromList(LinkedList *l, void * v){
  if(l->size == 0){
    return;
  }
  Node *n = l->head;
  Node *prev = NULL;
  while(n != NULL){
    if(n->value == v){
      if(prev == NULL){
        l->head = n->next;
      }
      else{
        prev->next = n->next;
      }
      free(n);
      l->size--;
      return;
    }
    prev = n;
    n = n->next;
  }
}

int containsVal(LinkedList *l, void * v){
  Node *n = l->head;
  while(n != NULL){
    if(n->value == v){
      return 1;
    }
    n = n->next;
  }
  return 0;
}

void destroyList(LinkedList *l){
  Node *n = l->head;
  while(n != NULL){
    Node *temp = n;
    n = n->next;
    free(temp);
  }
  l->head = NULL;
  l->tail = NULL;
  l->size = 0;
}