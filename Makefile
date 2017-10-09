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
CPPFLAGS += -I$(HOME)/codes/include
CPPFLAGS += -I/home/homer/work/scuff-em/src/libs/libscuff
LDFLAGS += -L$(LIBDIR) -Wl,-rpath,$(LIBDIR)
LDFLAGS += -L$(HOME)/codes/lib -L$(HOME)/lib
CXXFLAGS += -O3 -DHAVE_CONFIG_H
#CXXFLAGS += -ggdb -O0

HR_LIBS = -lscuff

RDL_LIBS=-lreadline -lncurses
HDF5_LIBS= -L$(HOME)/codes/hdf5-parallel/lib -lhdf5_hl -lhdf5 
LB_LIBS=-llapack -lopenblas -lgomp -lgfortran 
OTHER_LIBS = $(RDL_LIBS) $(HDF5_LIBS) $(LB_LIBS)

LIBS = $(HR_LIBS) $(OTHER_LIBS) -lpthread

##################################################
##################################################
##################################################
OBJS = scuff-spectrum.o 

tBeyn411:		tBeyn411.o ./libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

tlibBeyn:      		tlibBeyn.o ./libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 

scuff-spectrum:		$(OBJS) ./libBeyn.a
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 

tOmegaDerivatives:	tOmegaDerivatives.o
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

libBeyn.a:		libBeyn.o libBeyn.h
			ar r libBeyn.a libBeyn.o

clean:
		/bin/rm *.o *.a
