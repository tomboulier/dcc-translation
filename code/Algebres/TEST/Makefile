#
# test de SaprseColMat
#
# LD Nov 94

BIN= $(HOME)/bin
LIBPATH = -L$(HOME)/lib/ 
INCPATH = -I$(HOME)/include

CFLAGS = $(INCPATH) -O3
CPLUS = g++
CC = gcc

.SUFFIXES : .c .C .o

all : distpolytest InterpoTest test-sparse mul-test-sparse 

MYxgolden : MYxgolden.C MYcfunc.o
	$(CPLUS) MYxgolden.C MYcfunc.o -o $(BIN)/MYxgolden $(CFLAGS) $(LIBPATH)\
		-lsparsemat  -lm

xamoeba : xamoeba.C 
	$(CPLUS) xamoeba.C -o $(BIN)/xamoeba $(CFLAGS) $(LIBPATH)\
		-lsparsemat -lalp -lm

xgasdev : xgasdev.C 
	$(CPLUS) xgasdev.C -o $(BIN)/xgasdev $(CFLAGS) $(LIBPATH)\
		-lsparsemat -lalp -lm

distpolytest : distpolytest.C 
	$(CPLUS) distpolytest.C -o $(BIN)/distpolytest $(CFLAGS) $(LIBPATH)\
		-lsparsemat -lalp -lm

InterpoTest : InterpoTest.C 
	$(CPLUS) InterpoTest.C -o $(BIN)/InterpoTest $(CFLAGS) $(LIBPATH) -lm

test-sparse : test-sparse.C 
	$(CPLUS) test-sparse.C -o $(BIN)/test-sparse $(CFLAGS) $(LIBPATH)\
		-lsparsemat -lart3d -lalp -lm

mul-test-sparse : mul-test-sparse.C 
	$(CPLUS) mul-test-sparse.C -o $(BIN)/mul-test-sparse $(CFLAGS) $(LIBPATH)\
		-lsparsemat -lart3d -lalp -lm


clean : 
	\rm $(BIN)/test-sparse $(BIN)/mul-test-sparse *.o

.c.o:
	$(CC) $(CFLAGS) -c $<
