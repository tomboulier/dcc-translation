#
# Compilation of ALP v4
#
# EB 93 94 95 96
#

############################################################
#
#  L I B R A R Y - S P E C I F I C   P A R T
#
#

# CC = cxx
CC = g++

OBJS =	rvecbase.o permut.o coptions.o\
	rmatbase.o rmatops.o raffbase.o rsvecbase.o\
	basic.o ibox.o
NAME = alp

DATA = $(HOME)
INSTALLDIR = $(HOME)

DIR = .

INCPATH = -I$(DIR) -I$(DATA)/include  -I/usr/local/include

LIBPATH =  -L$(DATA)/lib -L.

CXXFLAGS = $(INCPATH) -O3

all : install

ref :
	latex2e alp-ref

INCS =	alp.h rvecbase.h rsvecbase.h rmatbase.h\
	rmatops.h raffbase.h permut.h coptions.h\
	basic.h ibox.h

############################################################
#
#  P R I V A T E   P A R T
#

LIB = lib$(NAME).a

# set to ar tvs for DEC Alpha ALPHA
# for sgi	ar cruvs $(LIB) $(OBJS)  
lib: $(OBJS)
	ar rc  $(LIB) $(OBJS)  

clean:
	/bin/rm -f *.o *~ $(LIB) *.dvi *.log *.aux

realclean: clean
	/bin/rm -f *.d

install: lib
	/bin/cp -p $(LIB) $(INSTALLDIR)/lib/
	if [ ! -d $(INSTALLDIR)/include/$(NAME) ]; then \
		/bin/mkdir $(INSTALLDIR)/include/$(NAME); \
	fi
	/bin/cp -p $(INCS) $(INSTALLDIR)/include/$(NAME)
	/bin/mv $(INSTALLDIR)/include/$(NAME)/$(NAME).h $(INSTALLDIR)/include

remove:
	/bin/rm -f $(INSTALLDIR)/lib/$(LIB)
	/bin/rm -f $(INSTALLDIR)/include/$(NAME).h
	/bin/rm -f $(INSTALLDIR)/include/$(NAME)/*.h
	-/bin/rmdir $(INSTALLDIR)/include/$(NAME)

############################################################
#
#  A U T O M A T I C   D E P E N D E N C I E S
#

.SUFFIXES: .d .cc

.cc.o: 
	$(CC) $(CXXFLAGS) -c $<

%.d : %.cc
	@echo "Building dependencies : $<"
	@sh -ec 'g++ -M $(CXXFLAGS) $< | sed '\''s/$*\.o/& $*\.d/g'\'' > $@'

#-include $(OBJS:.o=.d)
