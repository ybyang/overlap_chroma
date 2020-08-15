TOPDIR=/home/sunpeng1/soft/bin
CHROMA=/home/sunpeng1/soft/bin/chroma_gpu_double
CONFIG=./chroma-config

CXX=$(shell $(CONFIG) --cxx) 
CXXFLAGS=$(shell $(CONFIG) --cxxflags) 
#CXXFLAGS=-D_PDF_ $(shell $(CONFIG) --cxxflags) -I$(TOPDIR)/qla/include -I. $(MGCXXFLAGS) 
LDFLAGS=$(shell $(CONFIG) --ldflags) 
LIBS=$(shell $(CONFIG) --libs) 

HDRS=inline_eigen_maker.h \
        hwilson_eigenop.h \
        overlap_eigenop.h \
        quda_utils.h \
        inline_propagator_multi_eigen.h \
        chebyshev_coeff.h

OBJS=chroma.o \
     inline_eigen_maker.o \
     inline_propagator_multi_eigen.o

chroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -rf chroma  $(OBJS) *~
