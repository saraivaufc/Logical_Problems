#ifndef _SMV_STRING_H
#define _SMV_STRING_H
typedef struct string{
  struct assoc *link;
  char *text;
} string_rec,*string_ptr;

string_ptr find_string();
void init_string();
#endif