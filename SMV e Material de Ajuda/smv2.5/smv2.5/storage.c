#include <stdio.h>
#include <windows.h>
#include <sys/types.h>
#include "storage.h"
static char *addrlimit;
static char *addrfree;

/* this routine initializes the storage manager */
void init_storage()
{
#ifdef MACH
  mach_init();		/* needed to make sbrk() work */
#endif MACH
  /* addrfree points to the first free byte
     addrlimit points to the memory limit */
    addrfree = addrlimit = 0;//(char *) sbrk(0);
}

  static char* lastPtr =0;

/* get ALLOCSIZE more bytes from the O.S. */
static void getmore(int num)
{
  char *na;
  int count;
/*  fprintf(stderr,"Getting %d more bytes\n",ALLOCSIZE); */
  /*if(addrlimit != (char *)sbrk(0)){ /* in case someone else did sbrk */
   /* sbrk((4 - (sbrk(0) % 4)) % 4);
    addrfree = addrlimit = (char *)sbrk(0);
    if(((unsigned)addrlimit) % 4 != 0)
      rpterr("Failed to allocate %d bytes: addrlimit = %xH, na = %xH\n",
	     ALLOCSIZE,(int)addrlimit,(int)na);
  }
  if((na = (char *)sbrk(ALLOCSIZE)) != addrlimit)
    rpterr("Failed to allocate %d bytes: addrlimit = %xH, na = %xH\n",
	   ALLOCSIZE,(int)addrlimit,(int)na);*/  
  if (addrlimit-addrfree > 1024){
    GlobalReAlloc(lastPtr,addrfree-lastPtr,GPTR);
    addrstart -= addrfree-lastPtr;
  }
  //else
    //WastedMem += addrlimit-addrfree;//GPTR GMEM_FIXED

  count = (num > ALLOCSIZE)? num:ALLOCSIZE;
  if ((na = (char*)(GlobalAlloc(GPTR, count))) == NULL)//8byte align
	  rpterr("Failed to allocate %d bytes: not enough memory\n", ALLOCSIZE);
  addrfree = na;
  addrlimit = addrfree + count;
  lastPtr = na;  
  addrstart += count;
}

/* provide malloc for miscellaneuos storage allocation */
char *smv_malloc(n)
int n;
{
  if(n % 4)n=n+4-(n%4);  /* always allocate multiple of four bytes */
  if (addrfree + n > addrlimit)	  
	  getmore(n);
  
  {
    char *r = addrfree;
    addrfree += n;
    return(r);
  }
}

/* very simple implementation of free */
void smv_free(p)
char *p;
{
  return;
}

/* initialize a record manager.
   sets the free list to NULL,
   and makes sure the record size
   is at least big enough for a pointer */

mgr_ptr new_mgr(rec_size)
int rec_size;
{
  register mgr_ptr mp = (mgr_ptr)smv_malloc(sizeof(struct mgr));
  mp->free.link = 0;
  mp->rec_size = rec_size;
  mp->count = 0;
  mp->free_hook = 0;
  return(mp);
}



/* get a new record. if the free list
   is not empty, pull the first record off this
   list. else get ALLOCSIZE more bytes and
   make a new block of records. Link all
   of these record into a free list.
   then get the first element of this list
   by calling new_rec recursively. */

rec_ptr new_rec(mp)
register mgr_ptr mp;
{
  register rec_ptr p1;
  if(mp->free.link){
    rec_ptr r = mp->free.link;
    mp->free.link = (mp->free.link)->link;
    r->link = 0;
    return(r);
  }
  getmore(ALLOCSIZE);
  p1 = &(mp->free);
  while(addrlimit-addrfree >= mp->rec_size){
    p1->link = (rec_ptr)addrfree;
    p1 = (rec_ptr)addrfree;
    addrfree += mp->rec_size;
  }
  /* link field of the last record should be NIL (by Hiromi, 1998.5.5) */
  p1->link = ((rec_ptr)0);
  return(new_rec(mp));
}

/* put a record on the free list */
void free_rec(mp,r)
register mgr_ptr mp;
rec_ptr r;
{
  register rec_ptr rp = r;
  if(mp->free_hook)(*(mp->free_hook))(rp);
  rp->link = mp->free.link;
  mp->free.link = rp;
}
/* make a copy of r for mp, it is not in mp*/
rec_ptr dup_rec(mp,r)
mgr_ptr mp;
rec_ptr r;
{
  register rec_ptr res = new_rec(mp);
  //PC bcopy(r,res,mp->rec_size);
  memcpy(res, r,mp->rec_size);
  res->link = 0;
  return(res);
}

