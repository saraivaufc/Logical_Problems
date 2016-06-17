#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <sys/types.h>

#ifdef SMV_SIGNALS
#include <signal.h>
#endif

#include <time.h>
#include "storage.h"
#include "node.h"
#include "hash.h"
#include "bdd.h"
#include "assoc.h"
#include "y.tab.h"
#include <setjmp.h>
#include "lost.h"
#include <time.h>
#include <string.h>

#ifndef CLK_TCK
# define CLK_TCK 60
#endif

void init_eval(){}

int KEYTABLESIZE=16381, APPLY_CACHE_SIZE=16381, MINI_CACHE_SIZE=16381;

int heuristics = 0,verbose = 0,option_print_reachable=0,
  option_forward_search = 0,option_round_robin = 0,
  option_incremental = 0,option_restrict_trans = 0, option_interactive = 0,
  option_AG_only = 0, option_conj_part = 0;
int interactive_mode = 0;
int conj_part_limit = 0;
double fudge_factor;
char *output_order_file = 0,*input_order_file = 0;
static char *input_file = "stdin";
static char *myname = "smv";
static hash_ptr seq_hash;
//static char* addrstart;
int addrstart;//PC


extern node_ptr parse_tree;
#ifdef REORDER
int reorder = 0;
int reorder_bits = 10;
int reorder_size = 5000;
int reorder_maxsize = 300000;
#endif
#ifdef OTHER_SIMP
int option_othersimp = 0;
int option_gc_factor = 3;
int option_gc_limit = 1024; /* max size of memory in MB. Default is 1GB. */
int option_early = 0;
int option_checktrans = 0;
int option_drip=0;
int option_output_often=0;
#endif

#ifdef SMV_SIGNALS
/* Quit after outputting the order file with -quit option (default),
   but don't quit with -noquit. */
int option_quit = 1;

void signal_handler(sig)
int sig;
{
  int i;

  printf("\nSMV caught signal %d... ",sig);
  if(sig==10) /* SIGUSR1 */
    if(!reorder) {
      printf("Dynamic variable reordering ON at %d < nodes < %d\n\n",
	     reorder_size, reorder_maxsize);
      reorder = 1;
    }
    else {
      printf("Dynamic variable reordering OFF\n\n");
      reorder = 0;
    }
  else if(sig==12) { /* SIGUSR2 */
    printf("Force garbage collection\n\n");
    reset_maxnodes();
  }
  else {
    printf("\n");
    print_usage();
    exit(1);
  }
  for(i = 0 ; i < 32 ; i++)
    signal(i, signal_handler);
}
#endif

void main(argc,argv)
int argc;
char **argv;
{
  extern int yylineno;
#ifdef SMV_SIGNALS
  int i;

  for(i = 0 ; i < 32 ; i++)
    signal(i, signal_handler);

  { 
    FILE *ff = fopen(".smv-pid","w");
    fprintf(ff,"%d",getpid());
    fclose(ff);
  }
#endif
  init_storage();
  addrstart = 0;//addrstart = (char *) sbrk(0);
  argc--;
  myname = *(argv++);
  while(argc){
    if(argc == 1 && (**argv) != '-'){
      open_input(*(argv++));
      argc--;
      continue;
    }
    if(strcmp(*argv,"-rr")==0){
      argv++;
      argc--;
      option_round_robin = 1;
      continue;
    }
    if(strcmp(*argv,"-r")==0){
      argv++;
      argc--;
      option_print_reachable = 1;
      continue;
    }
    if(strcmp(*argv,"-f")==0){
      argv++;
      argc--;
      option_forward_search = 1;
      continue;
    }
    if(strcmp(*argv,"-inc")==0){
      argv++;
      argc--;
      option_incremental = 1;
      continue;
    }
    if(strcmp(*argv,"-int")==0){
      argv++;
      argc--;
      option_interactive = 1;
      //setlinebuf(stdout);
      continue;
    }
    if(strcmp(*argv,"-AG")==0){
      argv++;
      argc--;
      option_AG_only = 1;
      continue;
    }
    if(strcmp(*argv,"-h")==0){
      heuristics = 1;
      argv++;
      argc--;
      continue;
    }
#ifdef REORDER
    if(strcmp(*argv, "-reorder") == 0) {
      argv++;
      argc--;
      reorder = 1;
      continue;
    }
#endif
#ifdef OTHER_SIMP
    if(strcmp(*argv, "-early") == 0) {
      argv++;
      argc--;
      option_early = 1;
      continue;
    }
    if(strcmp(*argv, "-noearly") == 0) {
      argv++;
      argc--;
      option_early = 0;
      continue;
    }
    if(strcmp(*argv, "-checktrans") == 0) {
      argv++;
      argc--;
      option_checktrans = 1;
      continue;
    }
    if(strcmp(*argv, "-drip") == 0) {
      argv++;
      argc--;
      option_drip = 1;
      continue;
    }
#endif /* OTHER_SIMP */
#ifdef VERSION
    if(strcmp(*argv, "-version") == 0) {
      argv++;
      argc--;
      printf(VERSION);printf("\n");
      return(0);
    }
#endif
#ifdef SMV_SIGNALS
    if(strcmp(*argv, "-quit") == 0) {
      argv++;
      argc--;
      option_quit = 1;
      continue;
    }
    if(strcmp(*argv, "-noquit") == 0) {
      argv++;
      argc--;
      option_quit = 0;
      continue;
    }
#endif
    if(argc<2)rpterr("command line error");
    if(strcmp(*argv,"-v")==0){
      argv++;
      sscanf(*(argv++),"%d",&verbose);
      //setlinebuf(stdout);
    }
#ifdef OTHER_SIMP
    else if((strcmp(*argv, "-simp") == 0)
	    || (strcmp(*argv, "-othersimp") == 0)) {
      argv++;
      sscanf(*(argv++),"%d",&option_othersimp);
      if(option_othersimp < 0 || option_othersimp > 2) {
	fprintf(stderr,"Error: -othersimp %d: must be 0, 1 or 2\n");
	exit(1);
      }
    }
    else if(strcmp(*argv, "-gcfactor") == 0) {
      argv++;
      sscanf(*(argv++),"%d",&option_gc_factor);
    }
    else if(strcmp(*argv, "-gclimit") == 0) {
      argv++;
      sscanf(*(argv++),"%d",&option_gc_limit);
    }
    else if(strcmp(*argv,"-oo")==0){
      argv++;
      output_order_file = *(argv++);
      option_output_often = 1;
    }
#endif
    else if(strcmp(*argv,"-i")==0){
      argv++;
      input_order_file = *(argv++);
    }
    else if(strcmp(*argv,"-o")==0){
      argv++;
      output_order_file = *(argv++);
    }
    else if(strcmp(*argv,"-k")==0){
      argv++;
      sscanf(*(argv++),"%d",&KEYTABLESIZE);
    }
    else if(strcmp(*argv,"-c")==0){
      argv++;
      sscanf(*(argv++),"%d",&APPLY_CACHE_SIZE);
    }
    else if(strcmp(*argv,"-m")==0){
      argv++;
      sscanf(*(argv++),"%d",&MINI_CACHE_SIZE);
    }
    else if(strcmp(*argv,"-cp")==0){
      argv++;
      sscanf(*(argv++),"%d",&conj_part_limit);
      option_conj_part = 1;
    }
#ifdef REORDER
    else if(strcmp(*argv, "-reorderbits") == 0) {
      argv++;
      sscanf(*(argv++),"%d",&reorder_bits);
      reorder_bits *= 2;
    }
    else if(strcmp(*argv, "-reordersize") == 0) {
      argv++;
      sscanf(*(argv++),"%d",&reorder_size);
    }
    else if(strcmp(*argv, "-reordermaxsize") == 0) {
      argv++;
      sscanf(*(argv++),"%d",&reorder_maxsize);
    }
#endif
    else rpterr("undefined: %s",*argv);
    argc -= 2;
  }
#ifdef OTHER_SIMP
  /* Compute the actual amount of memory left for bdds */
  option_gc_limit -= (KEYTABLESIZE*sizeof(bdd_ptr)
		      + APPLY_CACHE_SIZE*sizeof(apply_rec))/(1024*1024);
#endif
  if(verbose){
    fprintf(stderr,"Key table size: %d\n",KEYTABLESIZE);
    fprintf(stderr,"Apply cache size: %d\n",APPLY_CACHE_SIZE);
    fprintf(stderr,"Variable ordering heuristics: ");
    if(heuristics)fprintf(stderr,"ON, factor = %g\n",fudge_factor);
    else fprintf(stderr,"OFF\n");
#ifdef OTHER_SIMP
    fprintf(stderr,"GC factor: %d, memory limit (adjusted): %dMB (%d BDD nodes)\n",
	    option_gc_factor,option_gc_limit,option_gc_limit*BDDS_IN_MB);
    if(option_early && option_AG_only && option_forward_search)
      fprintf(stderr,"Using early evaluation of AG specs\n");
    fprintf(stderr,"Simplification algorithm: ");
    if(option_othersimp==0)fprintf(stderr,"memory efficient (0)\n");
    else if(option_othersimp==1)fprintf(stderr,"time efficient (1)\n");
    else fprintf(stderr,"both (2) \n");
#endif
  }
  if(MINI_CACHE_SIZE > APPLY_CACHE_SIZE)
    rpterr("mini-cache-size (%d) is larger than the cache-size (%d)",
	   MINI_CACHE_SIZE, APPLY_CACHE_SIZE);
  init_string();
  init_assoc();
  init_node();
  init_bdd();
  init_eval();
  
  if(verbose){fprintf(stderr,"Parsing..."); fflush(stderr);}
  if(yyparse())my_exit(1);
  if(verbose)fprintf(stderr,"done.\n");
  
  close_input();

  build_symbols();
  my_exit(0);
}

void open_input(filename)
char *filename;
{
  extern int yylineno;
  extern FILE *yyin;
  input_file = filename;
  if(!(yyin = fopen(filename,"r")))
    rpterr("cannot open %s for input",filename);
  yylineno = 1;
}

void close_input()
{
  extern int yylineno;
  input_file = 0;
  yylineno = 0;
}

static node_ptr atom_stack=0;

void undefined(s)
node_ptr s;
{
  start_err();
  print_node(stderr,s);
  fprintf(stderr," undefined");
  finish_err();
}

void redefining(s)
node_ptr s;
{
  start_err();
  fprintf(stderr,"redefining ");
  print_node(stderr,s);
  finish_err();
}

void circular(s)
node_ptr s;
{
  start_err();
  fprintf(stderr,"recursively defined: ");
  print_node(stderr,s);
  finish_err();
}

void toomanyvars(s)
node_ptr s;
{
  start_err();
  fprintf(stderr,"too many variables");
  finish_err();
}

void start_err()
{
  extern int yylineno;
  fprintf(stderr,"\n");
  if(input_file)fprintf(stderr,"file %s: ",input_file);
  if(yylineno)fprintf(stderr,"line %d: ",yylineno);
}

jmp_buf longjmp_buf;
int longjmp_on_err = 0;
void finish_err()
{
  fprintf(stderr,"\n");
  while(atom_stack){
    node_ptr s = car(atom_stack);
    atom_stack = cdr(atom_stack);
    fprintf(stderr,"in definition of ");
    print_node(stderr,s);
    if(s->lineno)
      fprintf(stderr," at line %d",s->lineno);
    fprintf(stderr,"\n");
  }
  if(longjmp_on_err)longjmp(longjmp_buf,1);
  my_exit(1);
}

int my_setjmp()
{
  int v;
  longjmp_on_err = 1;
  v = setjmp(longjmp_buf);
  if(v)
    longjmp_on_err = 0;
  return(v);
}

void cancel_my_setjmp()
{
    longjmp_on_err = 0;
}  
  

void my_exit(n)
int n;
{
  if(verbose)fprintf(stderr,"%s: exit(%d)\n ",myname,n);
#ifdef SERGEYDEBUG
  if(verbose && n)print_usage();
#endif
  exit(n);
}

#ifdef SERGEYDEBUG
extern unsigned apply_cache_access_count;
extern unsigned apply_cache_hit_count;
#endif

void print_usage()
{
  //struct tms buffer;
  //struct tm buffer;
  clock_t finish;
  printf("\nresources used:\n");
  finish = clock();
  printf("processor time: %g s, \n", (double)(finish) / CLOCKS_PER_SEC);
  //times(&buffer);//PC
  //printf("user time: %g s, system time: %g s\n",
//	 buffer.tms_utime/(double)CLK_TCK,
//	 buffer.tms_stime/(double)CLK_TCK);
  printf("BDD nodes allocated: %d\n",get_bdd_nodes_allocated());
  printf("Bytes allocated: %d\n",addrstart);//PC (unsigned)((char *)sbrk(0)-addrstart));
#ifdef SERGEYDEBUG
  printf("Apply cache: %u, hits: %u, hit rate: %2.1f %%\n",
	 apply_cache_access_count, apply_cache_hit_count,
	 (100.0*(float)apply_cache_hit_count)/(float)apply_cache_access_count);
#endif
}

/*VARARGS1*/
rpterr(s,a1,a2,a3,a4)
char *s,*a1,*a2,*a3,*a4;
{
  start_err();
  fprintf(stderr,s,a1,a2,a3,a4);
  finish_err();
}

/*VARARGS1*/
void catastrophe(s,a1,a2,a3,a4)
char *s,*a1,*a2,*a3,*a4;
{
  fprintf(stderr,"\n\n*** internal error *** ");
  fprintf(stderr,s,a1,a2,a3,a4);
  fprintf(stderr,"\nPlease report this error to sergey.berezin@cs.cmu.edu (not to mcmillan@cs...)\n");
  fprintf(stderr,"Send a copy of this output and your input.\n");
  my_exit(1);
}

void push_atom(s)
node_ptr s;
{
  atom_stack = cons(s,atom_stack);
}

void pop_atom()
{
  node_ptr temp;
  if(!atom_stack)catastrophe("pop_atom: stack empty");
  temp = cdr(atom_stack);
  free_node(atom_stack);
  atom_stack = temp;
}

void yyerror(s)
char *s;
{
    extern yytext;
    start_err();
    fprintf(stderr,"at token \"%s\": %s\n",&yytext,s);
    if(!interactive_mode)finish_err();
}

yywrap()
{
  return(1);
}

int indent_size = 0;
void indent(stream)
FILE *stream;
{
  int i;
  for(i=0;i<indent_size;i++)fprintf(stream,"  ");
}

void indent_node(stream,s1,n,s2)
FILE *stream;
char *s1,*s2;
node_ptr n;
{
  indent(stream);
  fprintf(stream,"%s",s1);
  print_node(stream,n);
  fprintf(stream,"%s",s2);
}
