##################################################
##################################################
##################################################
SCUFFEM  = /home/homer/work/scuff-em-installation
#SCUFFEM  = /home/homer/work/scuff-em-debug
LIBDIR   = $(SCUFFEM)/lib

##################################################
##################################################
##################################################
CPPFLAGS += -I$(SCUFFEM)/include/scuff-em -I$(HOME)/include -I. -fopenmp
CPPFLAGS += -I$(HOME)/codes/include -I$(HOME)/work/libBeyn
CPPFLAGS += -I/home/homer/work/scuff-em/src/libs/libscuff
LDFLAGS += -L$(LIBDIR) -Wl,-rpath,$(LIBDIR)
LDFLAGS += -L$(HOME)/codes/lib -L$(HOME)/lib -L$(HOME)/work/libBeyn
CXXFLAGS += -O3 -DHAVE_CONFIG_H
#CXXFLAGS += -ggdb -O0

CPPFLAGS += -I/home/home/work/libBeyn

HR_LIBS = -lBeyn -lscuff

RDL_LIBS=-lreadline -lncurses
HDF5_LIBS= -lhdf5_hl -lhdf5 
LB_LIBS=-llapack -lopenblas -lgomp -lgfortran 
OTHER_LIBS = $(RDL_LIBS) $(HDF5_LIBS) $(LB_LIBS)

LIBS = $(HR_LIBS) $(OTHER_LIBS) -lpthread

##################################################
##################################################
##################################################
OBJS = scuff-spectrum.o 

scuff-spectrum:		$(OBJS) libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 

tOmegaDerivatives:	tOmegaDerivatives.o
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

libBeyn.a:		libBeyn.o libBeyn.h
			ar r libBeyn.a libBeyn.o

clean:
		/bin/rm *.o *.a
