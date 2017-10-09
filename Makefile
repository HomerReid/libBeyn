##################################################
# set path to SCUFF-EM installation
##################################################
SCUFFEM  = /home/homer/work/scuff-em-installation

SCUFFINC  =  $(SCUFFEM)/include/scuff-em
SCUFFLIB  =  $(SCUFFEM)/lib
CPPFLAGS  += -I$(SCUFFINC)
LDFLAGS   += -L$(SCUFFLIB) -Wl,-rpath,$(SCUFFLIB)
SCUFF_LIBS = -lscuff

##################################################
##################################################
##################################################
CPPFLAGS+=-I$(HOME)/include -I$(HOME)/codes/include -I. -fopenmp
LDFLAGS+=-L$(HOME)/codes/lib -L$(HOME)/lib

CXXFLAGS+=-O3 -DHAVE_CONFIG_H
#CXXFLAGS+=-ggdb -O0

HDF5_LIBS  =  -L$(HOME)/codes/hdf5-parallel/lib -lhdf5_hl -lhdf5 
LB_LIBS    = -llapack -lopenblas -lgomp -lgfortran 
OTHER_LIBS = $(RDL_LIBS) $(HDF5_LIBS) $(LB_LIBS)

LIBS = $(SCUFF_LIBS) $(OTHER_LIBS) -lpthread

##################################################
##################################################
##################################################
all:			tBeyn411 scuff-spectrum

tBeyn411:		tBeyn411.o ./libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

scuff-spectrum:		scuff-spectrum.o ./libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 

libBeyn.a:		libBeyn.o libBeyn.h
			ar r libBeyn.a libBeyn.o

clean:
		/bin/rm *.o *.a
