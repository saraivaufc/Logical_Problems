#ifndef _NODE_H
#define _NODE_H

typedef union {
  int inttype;
  struct node *nodetype;
  struct string *strtype;
  struct bdd *bddtype;
} value;

typedef struct node{
  struct node *link;
  short int type,lineno;
  value left,right;
} node_rec,*node_ptr;

node_ptr new_node(),find_node();
void init_node(),free_node();
void print_node();
#ifdef OTHER_SIMP
void free_list();
#endif
node_ptr subst_node(),map(),key_node();
node_ptr cons(),car(),cdr(),append(),reverse();
node_ptr unify_node();

#define NIL ((node_ptr)0)
#define FAILURE_NODE ((node_ptr)(-1))
#endif