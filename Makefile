TOPDIR=/public/home/sunpeng/soft/bin
CHROMA=/public/home/sunpeng/soft/bin/chroma_gpu_double
CONFIG=./chroma-config

CXX=$(shell $(CONFIG) --cxx) 
CXXFLAGS=$(shell $(CONFIG) --cxxflags) "-w"
LDFLAGS=$(shell $(CONFIG) --ldflags) 
LIBS=$(shell $(CONFIG) --libs) 

HDRS=inline_eigen_maker.h \
        hwilson_eigenop.h \
        overlap_eigenop.h \
        quda_utils.h \
        inline_propagator_multi_eigen.h \
        chebyshev_coeff.h \
		readwrite.h

OBJS=chroma.o \
     inline_eigen_maker.o \
     inline_propagator_multi_eigen.o \
	 readwrite.o

chroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

BUILD_TYPE=build_quda
PWD_DIR=$(shell pwd)

build:
	cd $(BUILD_TYPE) && TMPDIR=$(PWD_DIR)/$(BUILD_TYPE)/tmp make -j 8 VERBOSE=1 2>&1 | tee build.log

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -rf chroma  $(OBJS) *~
