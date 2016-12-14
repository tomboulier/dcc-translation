
# installation sur MedEtPhy 
# mai 2005
# ########################## #### #### ####
LIBROOT = $(HOME)
LIBPATH = -L$(HOME)/lib  
INCPATH = -I$(HOME)/include 

CFLAGS = -O3 $(INCPATH) 
CPLUS = g++
CC = $(CPLUS) $(CFLAGS)

BINDIR = $(LIBROOT)/bin


.SUFFIXES : .c .C .o

all : GFBLineSimulData GFBCircularSimulData

GFBLineSimulData : GFBLineSimulData.C 
	$(CPLUS) GFBLineSimulData.C -o $(BINDIR)/GFBLineSimulData $(CFLAGS) $(LIBPATH)\
		-ltomographie  -lrealimage -lm 

GFBCircularSimulData : GFBCircularSimulData.C 
	$(CPLUS) -w GFBCircularSimulData.C -o $(BINDIR)/GFBCircularSimulData $(CFLAGS) $(LIBPATH)\
		-ltomographie  -lrealimage -lsparsemat -lm 

clean : 
	/bin/rm *.o *~ $(BINDIR)/FBfbi $(BINDIR)/fbi

.C.o:
	$(CPLUS) $(CFLAGS) -c $<

