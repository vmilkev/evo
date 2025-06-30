APPNAME := evolm

SRCDIR := src
BDIR_PYLIB := release
BDIR_CLIB := release2
BDIR_APP := release3
TESTDIR := tests
BDIR_TST := $(TESTDIR)/release
TSRCDIR := $(TESTDIR)/src

EXTSRC2 := ../libplinkio/src
EXTLIBNAME := plinkio

# MAC
LIBRARY_PATH := /Users/au383883/Documents/MY/codebase/evo/libplinkio/release_lib
LIBRARY_PATH_2 := /Users/au383883/Documents/MY/codebase/evo/evolm/release2
# DEBIAN
#LIBRARY_PATH := /usr/home/qgg/vimi/codebase/evo/libplinkio/release_lib
#LIBRARY_PATH_2 := /usr/home/qgg/vimi/codebase/evo/evolm/release2
# UBUNTU
#LIBRARY_PATH := /home/vimi/Documents/codebase/evo/libplinkio/release_lib
#LIBRARY_PATH_2 := /home/vimi/Documents/codebase/evo/evolm/release2

# determine OS and source of BLAS lib
UNAME_OS := $(shell uname -s)
SOURSE_BLAS_IS_INTEL := 0

# determine compiler
ifeq ($(UNAME_OS),Darwin)
	CXX := clang++
	LIBRARY_PATH := /Users/au383883/Documents/MY/codebase/evo/libplinkio/release_lib
	LIBRARY_PATH_2 := /Users/au383883/Documents/MY/codebase/evo/evolm/release2
else
# cluster
#LIBRARY_PATH := /usr/home/qgg/vimi/codebase/evo/libplinkio/release_lib
#LIBRARY_PATH_2 := /usr/home/qgg/vimi/codebase/evo/evolm/release2
# Ubuntu
	LIBRARY_PATH := /home/vimi/Documents/codebase/evo/libplinkio/release_lib
	LIBRARY_PATH_2 := /home/vimi/Documents/codebase/evo/evolm/release2
	SHELL := /bin/bash
	export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):$(LIBRARY_PATH):$(LIBRARY_PATH_2)
	
	ifdef SOURSE_BLAS_IS_INTEL
		CXX := icpx
	else
		CXX := g++
	endif
endif

PWDADDR := $(shell pwd)
PYBINDINCL := $(shell python3 -m pybind11 --includes)
PYSUFF := $(shell python3-config --extension-suffix)
VERSUFF := .so

ifeq ($(UNAME_OS),Darwin) # using apple's Accelerate
	CFLAGS_TEST :=  -I$(SRCDIR) $(PYBINDINCL) -I$(EXTSRC2)
	CFLAGS_APP := -I$(SRCDIR) -I$(EXTSRC2)
	CFLAGS_PYLIB := -DPYBIND $(PYBINDINCL) -I$(EXTSRC2)
	CFLAGS_CLIB := -I$(EXTSRC2)

	EXTRA_CFLAGS := -Wall -Wextra -pedantic -Wno-c99-extensions -D UTEST -O3 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/libomp/include -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++17 -pthread

	LIBFLAGS := -undefined dynamic_lookup -shared -lpthread -lm -lomp -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/openblas/lib -lopenblas
	LIBFLAGS_TEST := -undefined dynamic_lookup -lpthread -lm -lomp -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/openblas/lib -lopenblas
	EXTRA_LIBFLAGS := -L$(LIBRARY_PATH) -l$(EXTLIBNAME)
	EXTRA_LIBFLAGS_TEST := -L$(BDIR_CLIB) -l$(APPNAME)
else
	ifeq ($(UNAME_OS),Linux)
#		ifeq ( $(SOURSE_BLAS_IS_INTEL), 1) # using Intel's MKL
		ifdef SOURSE_BLAS_IS_INTEL
			CFLAGS_TEST := -O3 -qopenmp -I$(SRCDIR) -I$(EXTSRC2) -DMKL_ILP64 -I${MKLROOT}/include $(PYBINDINCL) -std=c++17
			CFLAGS_APP := -O3 -qopenmp -I$(SRCDIR) -I$(EXTSRC2) -DMKL_ILP64 -I${MKLROOT}/include -std=c++17
			CFLAGS_PYLIB := -O3 -qopenmp -I$(EXTSRC2) -DPYBIND -DMKL_ILP64 -I${MKLROOT}/include $(PYBINDINCL) -fPIC -std=c++17
			CFLAGS_CLIB := -O3 -qopenmp -I$(EXTSRC2) -DMKL_ILP64 -I${MKLROOT}/include -fPIC -std=c++17
			
			EXTRA_CFLAGS := -D UTEST -D intelmkl -Wall -pthread
			EXTRA_LIBFLAGS := -L$(LIBRARY_PATH) -l$(EXTLIBNAME)
			EXTRA_LIBFLAGS_TEST := -L$(BDIR_CLIB) -l$(APPNAME)
			EXTRA_LIBFLAGS_2 :=
			
			LIBFLAGS := -shared -qopenmp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
			LIBFLAGS_TEST := -qopenmp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
		else # using openBLAS
			CFLAGS_APP := -I$(SRCDIR) -I(EXTSRC) -I$(EXTSRC2)
			CFLAGS_PYLIB := -DPYBIND $(PYBINDINCL) -I(EXTSRC) -I$(EXTSRC2) -fPIC
			EXTRA_CFLAGS := -Wall -Wextra -pedantic -D UTEST -03 -fopenmp -I/opt/openBLAS/include -std=c++17
			LIBFLAGS_PY := -shared -lpthread -lm -lgomp /opt/openBLAS/lib/libopenblas.a
			LIBFLAGS_TEST := -lpthread -lm -lgomp /opt/openBLAS/lib/libopenblas.a

			EXTRA_LIBFLAGS := -L$(LIBRARY_PATH) -l$(EXTLIBNAME)
			EXTRA_LIBFLAGS_TEST := -L$(BDIR_CLIB) -l$(APPNAME)
		endif
	endif
endif

ifdef perf_num
	EXTRA_CFLAGS += -D PERF_TEST
endif

TARGET_PYLIB := $(BDIR_PYLIB)/$(APPNAME)$(PYSUFF)
TARGET_CLIB := $(BDIR_CLIB)/lib$(APPNAME)$(VERSUFF)
TARGET_APP := $(BDIR_APP)/$(APPNAME)_app.exe
TARGET_TEST := $(BDIR_TST)/$(APPNAME)_test.exe

CXX: $(TARGET_PYLIB) $(TARGET_CLIB) $(TARGET_TEST) $(TARGET_APP)

SOURCES_ALL := $(wildcard $(SRCDIR)/*.cpp)
SOURCES = $(filter-out $(SRCDIR)/evo_profiling.cpp, $(wildcard $(SRCDIR)/*.cpp))
OBJECT_PYL := $(patsubst $(SRCDIR)/%,$(BDIR_PYLIB)/%,$(SOURCES:.cpp=.o))
OBJECT_LIB := $(filter-out $(BDIR_CLIB)/wrapper.o, $(patsubst $(SRCDIR)/%,$(BDIR_CLIB)/%,$(SOURCES:.cpp=.o)))
OBJECT_APP := $(filter-out $(BDIR_APP)/wrapper.o, $(patsubst $(SRCDIR)/%,$(BDIR_APP)/%,$(SOURCES_ALL:.cpp=.o)))

SOURCES_TST := $(wildcard $(TESTDIR)/$(SRCDIR)/*.cpp)
OBJECT_TST := $(patsubst $(TSRCDIR)/%,$(BDIR_TST)/%,$(SOURCES_TST:.cpp=.o))

#----------------------------------------------------------------

# build lib for python
$(TARGET_PYLIB): $(OBJECT_PYL)
	@echo " Linking python lib ..."
	$(CXX) $^ $(LIBFLAGS) $(EXTRA_LIBFLAGS) -o $(TARGET_PYLIB)

# build catch2 test
$(TARGET_TEST): $(OBJECT_TST)
	@echo " Linking tests ..."
	$(CXX) $^ $(LIBFLAGS_TEST) $(EXTRA_LIBFLAGS) $(EXTRA_LIBFLAGS_TEST) -o $(TARGET_TEST)

# build shared c++ lib
$(TARGET_CLIB): $(OBJECT_LIB)
	@echo " Linking shared lib ..."
	$(CXX) $^ $(LIBFLAGS) $(EXTRA_LIBFLAGS) -o $(TARGET_CLIB)

# build c++ app
$(TARGET_APP): $(OBJECT_APP)
	@echo " Linking c++ app ..."
	$(CXX) $^ $(LIBFLAGS_TEST) $(EXTRA_LIBFLAGS) -o $(TARGET_APP)

#----------------------------------------------------------------

# Compiling python lib
$(BDIR_PYLIB)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling pylib ..."
	@mkdir -p $(BDIR_PYLIB)
	$(CXX) -c $(CFLAGS_PYLIB) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling c++ lib
$(BDIR_CLIB)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling clib ..."
	@mkdir -p $(BDIR_CLIB)
	$(CXX) -c $(CFLAGS_CLIB) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling catch2 tests
$(BDIR_TST)/%.o: $(TSRCDIR)/%.cpp
	@echo " Compilling catch2 tests ..."
	@mkdir -p $(BDIR_TST)
	$(CXX) -c $(CFLAGS_APP) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling c++ app
$(BDIR_APP)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling c++ app ..."
	@mkdir -p $(BDIR_APP)
	$(CXX) -c $(CFLAGS_APP) $(EXTRA_CFLAGS) -o $@ $< 

#----------------------------------------------------------------

clean:
	@echo " Cleaning..."
	@rm -fr $(BDIR_PYLIB) $(BDIR_CLIB) $(BDIR_TST) $(BDIR_APP) $(TARGET_PYLIB) $(TARGET_TEST) 2>/dev/null || true
	@rm amatrix_* py_out_solution_model_* cpp_solution_model_* 2>/dev/null || true

rclean:
	@echo " Cleaning..."
	@rm amatrix_* py_out_solution_model_* cpp_solution_model_* model_matrix_* sparse_solution_* solution_* *.log *.dmbin *.corbin 2>/dev/null || true

test:
	@echo " Integrated tests with Python interface ..."
	$(TESTDIR)/$(SRCDIR)/launch_test_models.sh
ifeq ($(UNAME_OS),Darwin)
	@echo " Unit tests ..."
	@export DYLD_LIBRARY_PATH=$(LIBRARY_PATH); $(TARGET_TEST)
else
	$(TARGET_TEST)
endif

run:
	@echo " Running ..."
ifeq ($(UNAME_OS),Darwin)
	@export DYLD_LIBRARY_PATH=$(LIBRARY_PATH); time $(TARGET_APP)
else
	time $(TARGET_APP)
endif

cleanlogs:
	@echo " Cleaning logs ..."
	@rm *.log

.PHONY: clean
