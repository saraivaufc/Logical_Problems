#ifndef _LOST_H
#define _LOST_H
#include <stdio.h>
#include "node.h"
#include "smvstring.h"
#include "bdd.h"

int yylook();
int yywrap();
string_ptr find_string(char* x);
void indent_node(FILE * stream,char *s1,node_ptr n,char *s2);
//void rpterr();//char *s,char *a1,char*a2,char*a3,char*a4);
int yyback();

void catastrophe();
void type_error();
void bzero(void *s, int n);
void set_variable_names();
void output_order();
int yylex();
void yyerror(char *s);
node_ptr find_atom(node_ptr a);
void check_spec(node_ptr the_spec);
void compute_bound(node_ptr the_spec);
void goto_state(node_ptr s);
void assign_command(node_ptr var, node_ptr val);
void single_step();
void eval_command(node_ptr exp);
void trans_command(node_ptr n);
void init_command(node_ptr n);
void fair_command(node_ptr n);
void reset_command();
void build_symbols();
static void swapwords(bdd_ptr* a,bdd_ptr* b);
static void markbddlist(struct node *bddlist);
void check_bdd(bdd_ptr d);
void pr_status();
void yyoutput(int c);
void yyunput(int c);
void open_input(char *filename);
void init_assoc();
int yyparse();
void my_exit(int n);
void cancel_my_setjmp();
void close_input();
void undefined(node_ptr s);
void redefining(node_ptr s);
void circular(node_ptr s);
void toomanyvars();
void start_err();
void finish_err();
int get_bdd_nodes_allocated();
void print_usage();
void push_atom(node_ptr s);
void pop_atom();
int sprint_node(char *str,int size, node_ptr n);
int print_node_atcol(FILE *stream,node_ptr n,int col);
static void notanumber(node_ptr n);
void indent(FILE *stream);
void mygarbage();
static void range_error(node_ptr n);
void walk_leaves(void (*f)(),bdd_ptr d);
int value_bdd(bdd_ptr a);
bdd_ptr if_then_bdd(bdd_ptr a,bdd_ptr b);
bdd_ptr varset_diff(bdd_ptr a, bdd_ptr b);
void walk(void (*f)(),node_ptr l);
bdd_ptr collapse_vars(bdd_ptr a,bdd_ptr b,bdd_ptr v);
bdd_ptr collapse_no_shift(bdd_ptr a,bdd_ptr b);
void restart_bdd();
static bdd_ptr aux_make_var_mask(bdd_ptr d,int n,int l);
void reorder_variables();
bdd_ptr make_var_mask(bdd_ptr d);
void reset_maxnodes();
void report_and_exit();

























 
 




#endif