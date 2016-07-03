VERSION = 2.5.0
VERFLAG = -DVERSION="\"SMV $(VERSION). Sat May 9, 1998.\""

# build SMV with the dynamic variable reordering option
# comment out next line if not desired
REORDER = -DREORDER

# build SMV with the option of computing minimum and maximum path lengths
# comment out next two lines if not desired
TIMING = -DTIMING
TM4 = timing.m4

# build SMV on the MACH operating system
# comment out next line otherwise
# MACHLIB = -lmach

OPTFLAGS = -O2 
DFLAGS = -DBETTERHASH -DSMV_SIGNALS -DOTHER_SIMP
CFLAGS =  $(OPTFLAGS) -I. $(DFLAGS) $(REORDER) $(TIMING) $(VERFLAG)
LIB = -lm $(MACHLIB)
OBJ = grammar.o input.o main.o assoc.o bdd.o  hash.o node.o storage.o \
	string.o symbols.o
EXEC = smv
CC = gcc

# Pick the one that works under your system
# M4 = m4
M4 = gm4

YACC = yacc
LEX = lex

# On some Linux platforms you may have a new version of lex that
# requires a command line option for backward compatibility:
LEXFLAGS = -l
#LEXFLAGS = 

$(EXEC): $(OBJ) $(DEBUG)
	rm -f $(EXEC)
	$(CC) $(CFLAGS) $(OBJ) $(LIB) $(DEBUG) -o $(EXEC)
grammar.y:	grammar.m4
	rm -f grammar.y
	$(M4) $(TM4) grammar.m4 > grammar.y
grammar.c y.tab.h: grammar.y
	$(YACC) -d grammar.y; mv y.tab.c grammar.c
input.lex:	input.m4
	rm -f input.lex
	$(M4) $(TM4) input.m4 > input.lex
input.c: input.lex
	$(LEX) $(LEXFLAGS) input.lex; mv lex.yy.c input.c

main.o: main.c y.tab.h assoc.h  storage.h bdd.h hash.h node.h string.h
assoc.o: assoc.c assoc.h storage.h hash.h node.h
bdd.o: bdd.c assoc.h storage.h bdd.h hash.h node.h
symbols.o : symbols.c storage.h node.h hash.h bdd.h assoc.h y.tab.h
hash.o: hash.c storage.h hash.h
node.o: node.c y.tab.h assoc.h storage.h bdd.h hash.h node.h
storage.o: storage.c storage.h
string.o: string.c assoc.h storage.h hash.h string.h
grammar.o: grammar.c storage.h node.h hash.h assoc.h
input.o: input.c node.h y.tab.h

clean:
	rm -f *.o grammar.c grammar.y input.lex input.c y.tab.h
SRC= \
  makefile \
  README \
  assoc.c \
  assoc.h \
  bdd.c \
  bdd.h \
  examples/counter.smv \
  examples/dme1.smv \
  examples/dme2.smv \
  examples/gigamax.smv \
  examples/ring.smv \
  examples/semaphore.smv \
  examples/short.smv \
  examples/mutex.smv \
  examples/mutex1.smv \
  examples/periodic.smv \
  examples/syncarb5.smv \
  grammar.m4 \
  hash.c \
  hash.h \
  input.m4 \
  main.c \
  node.c \
  node.h \
  smv.1 \
  storage.c \
  storage.h \
  string.c \
  string.h \
  symbols.c \
  timing.m4 \
  NEW

tar:
	@srcs=""; \
	dir=`pwd`; dir=`basename $$dir`; \
	for f in $(SRC); \
	do \
	  srcs="$$srcs $${dir}/$${f}"; \
	done; \
	echo $$srcs; \
	cd ..; \
	tarfile="$${dir}.r$(VERSION).tar.gz"; \
	tar cf - $${srcs} | gzip -c > $$tarfile; \
	echo "> $$tarfile"; \
	mv $$tarfile $${dir}/; cd $${dir}
