##################################################
# set path to SCUFF-EM installation
##################################################
#SCUFFEM  = /usr/local
SCUFFEM  = ${SEI}/..

SCUFFINC    =  $(SCUFFEM)/include/scuff-em
SCUFFLIB    =  $(SCUFFEM)/lib
SCUFFLIB64  =  $(SCUFFEM)/lib64
CPPFLAGS   += -I$(SCUFFINC)
LDFLAGS    += -L$(SCUFFLIB) -Wl,-rpath,$(SCUFFLIB)
LDFLAGS    += -L$(SCUFFLIB64) -Wl,-rpath,$(SCUFFLIB64)
SCUFF_LIBS = -lscuff

##################################################
##################################################
##################################################
CPPFLAGS+=-I$(HOME)/include -I$(HOME)/codes/include -I. -fopenmp
LDFLAGS+=-L$(HOME)/codes/lib -L$(HOME)/lib

CXXFLAGS+=-O3

LB_LIBS    = -llapack -lopenblas -lgomp -lgfortran
OTHER_LIBS = $(RDL_LIBS) $(LB_LIBS)

LIBS = $(SCUFF_LIBS) $(OTHER_LIBS) -lpthread

##################################################
##################################################
##################################################
all:			tBeyn411

tBeyn411:		tBeyn411.o ./libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

libBeyn.a:		libBeyn.o libBeyn.h
			ar r libBeyn.a libBeyn.o

clean:
		/bin/rm *.o *.a
