#ifndef _BDD_H
#define _BDD_H

typedef struct bdd {
   unsigned dfield;
   struct bdd *left, *right, *next;
} bdd_rec, *bdd_ptr;

#ifndef NULL
#define NULL 0
#endif

extern bdd_ptr ZERO,ONE;

/* Data structure for the key table */
typedef struct {
   int n;
   int elements_in_table;
   bdd_ptr *hash_table_buf;
} keytable_rec, *keytable_ptr;

/* Data structure for apply record */
typedef struct {
  int op;
  bdd_ptr arg1,arg2,res;
} apply_rec;



/* Bits packed as follows:                                   */
/* 19:0  --- id    (0 to IDMASK : depth first id)            */
/* 30:20 --- level (0 to LEAFLEVEL)                          */
/* 31    --- mark  (TRUE or FALSE : used in graph traversal) */

#define IDMASK    0X0000FFFF /* This is never used. Reduced */
#define LEVELMASK 0X7FFF0000 /* Increased for 4 bits to 32K levels */
#define MARKMASK  0X80000000

#define IDLOW    0
#define LEVELLOW 16  /* This has changed as well */
#define MARKLOW  31

#define LEAFLEVEL 0X7FFF /* ... and this. */
#define ISLEAF(d) (GETLEVEL(d) == LEAFLEVEL)

#define GETFIELD(var, mask, low) ((int) ((var & mask) >> low))
#define PUTFIELD(var, val, mask, low) \
        (var = (var & ~(mask)) | (((unsigned) val) << low))

     /* #define GETID(d)        ((d)->dfield & IDMASK)
      * #define SETID(d, idval) (PUTFIELD((d)->dfield, idval, IDMASK, IDLOW))
      */

#define GETLEVEL(d) (GETFIELD((d)->dfield, LEVELMASK, LEVELLOW))
#define SETLEVEL(d, lval) (PUTFIELD((d)->dfield, lval, LEVELMASK, LEVELLOW))

#define SETMARK(d)   ((d)->dfield |= MARKMASK)
#define CLEARMARK(d) ((d)->dfield &= ~MARKMASK)
#define TESTMARK(d)  (((d)->dfield & MARKMASK) != 0X0)


#define AND_OP 1
#define OR_OP 2
#define XOR_OP 3
#define FORSOME_OP 4
#define NEXT_OP 11
#define PREV_OP 12
#define COMP_OP 6
#define SIMP_OP 7
#define COUNT_OP 8
#define ELIM_OP 9
#define SATISFY_OP 10

#ifdef OTHER_SIMP

#define SIMP_OP2 13
#define USE_BIG_CACHE 14
#define BDDS_IN_MB (1024*1024/sizeof(bdd_rec))

#else /* not OTHER_SIMP */

#define USE_BIG_CACHE 11

#endif /* OTHER_SIMP */

/* functions we provide */

void init_bdd();
bdd_ptr find_bdd();
void sweep_reduce();
void save_apply();
void insert_apply();
bdd_ptr find_apply();
void flush_apply();
void repairmark();
void renumber();
int size_bdd();
void mark_bdd();
bdd_ptr and_bdd();
bdd_ptr or_bdd();
bdd_ptr xor_bdd();
bdd_ptr not_bdd();
bdd_ptr forsome();
#ifdef OTHER_SIMP
bdd_ptr simplify_assuming2();
#endif
bdd_ptr simplify_assuming();
bdd_ptr sat_bdd();
double count_bdd(),n_count_bdd();
bdd_ptr save_bdd();
void release_bdd();
bdd_ptr leaf_bdd(),atomic_bdd(),r_shift(),f_shift(),r_collapse(),collapse();
bdd_ptr apply_bdd(),if_then_else_bdd();

#define IS_CURRENT_VAR(s) (((s)&1)==0)
#define IS_NEXT_VAR(s) (((s)&1)==1)
#define THE_CURRENT_VAR(s) (((s)<<1))
#define THE_NEXT_VAR(s) (((s)<<1)+1)
#define VAR_NUM(s) ((s)>>1)
#define NEXT_TO_CURRENT(s) ((s)-1)
#define CURRENT_TO_NEXT(s) ((s)+1)
#endif