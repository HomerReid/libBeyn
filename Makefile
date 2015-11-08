##################################################
##################################################
##################################################
SCUFFEM  = /home/homer/work/scuff-em-installation
#SCUFFEM  = /home/homer/work/scuff-em-debug
LIBDIR   = $(SCUFFEM)/lib
#LIBDIR = /home/homer/work/scuff-em-debug/src/libs/libscuff/.libs

##################################################
##################################################
##################################################
CPPFLAGS += -I$(SCUFFEM)/include/scuff-em -I$(HOME)/include -I. -fopenmp
CPPFLAGS += -I$(HOME)/codes/include
CPPFLAGS += -I/home/homer/work/scuff-em/src/libs/libscuff
LDFLAGS += -L$(LIBDIR) -Wl,-rpath,$(LIBDIR)
LDFLAGS += -L$(HOME)/codes/lib -L$(HOME)/lib 
CXXFLAGS += -O3 -DHAVE_CONFIG_H
#CXXFLAGS = -ggdb -O0

HR_LIBS=-lscuff

##################################################
##################################################
##################################################
CPPFLAGS += -I$(HOME)/work/libDCUTRI
LDFLAGS += -L$(HOME)/work/libDCUTRI
HR_LIBS=-lscuff -lDCUTRI

RDL_LIBS=-lreadline -lncurses
HDF5_LIBS= -lhdf5_hl -lhdf5 
LB_LIBS=-llapack -lopenblas -lgomp -lgfortran 
OTHER_LIBS = $(RDL_LIBS) $(HDF5_LIBS) $(LB_LIBS)

LIBS = $(HR_LIBS) -lnlopt $(OTHER_LIBS)

##################################################
##################################################
##################################################
OBJS = scuff-spectrum.o ObjectiveFunction.o

scuff-spectrum:		$(OBJS)
			$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		/bin/rm *.o *.a
