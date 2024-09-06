APPNAME := evoped

# set the project directories
SRCDIR := src
BDIR_PYLIB := release
BDIR_APP := release2
TESTDIR := tests
BDIR_TST := $(TESTDIR)/release
TSRCDIR := $(TESTDIR)/src
EXTSRC := ../evolm/src

EXTSRC2 := ../libplinkio/src
EXTLIB := ../libplinkio/release_lib
EXTLIBNAME := plinkio

# Mac
LIBRARY_PATH := /Users/au383883/Documents/MY/codebase/evo/libplinkio/release_lib
# Debian
#LIBRARY_PATH := /usr/home/qgg/vimi/codebase/evo/libplinkio/release_lib
# Ubuntu
#LIBRARY_PATH := /home/vimi/Documents/codebase/evo/libplinkio/release_lib

# determine OS and source of BLAS lib
UNAME_OS := $(shell uname -s)
SOURSE_BLAS_IS_INTEL := 1

# determine compiler
ifeq ($(UNAME_OS),Darwin)
	CXX := clang++
else
	SHELL := /bin/bash
	export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):$(LIBRARY_PATH)
	ifdef SOURSE_BLAS_IS_INTEL
		CXX := icpx
	else
		CXX := g++
	endif
endif

PWDADDR := $(shell pwd)
PYBINDINCL := $(shell python3 -m pybind11 --includes)
PYSUFF := $(shell python3-config --extension-suffix)

#$(info $$UNAME_OS is [${UNAME_OS}])
#$(info $$CXX is [${CXX}])
$(info $$UNAME_OS is [${UNAME_OS}])
#$(info $$SOURSE_BLAS_IS_INTEL is [${SOURSE_BLAS_IS_INTEL}])

ifeq ($(UNAME_OS),Darwin)
$(info compiling for OpenBlas libs on MacOS)
	CFLAGS_APP := -I$(SRCDIR) -I$(EXTSRC) -I$(EXTSRC2)
	CFLAGS_PYLIB := -DPYBIND $(PYBINDINCL) -I$(EXTSRC) -I$(EXTSRC2)

	EXTRA_CFLAGS := -Wall -Wextra -pedantic -Wno-c99-extensions -D UTEST -O3 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/libomp/include -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++17 -pthread

	LIBFLAGS_PY := -undefined dynamic_lookup -shared -lpthread -lm -lomp -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/openblas/lib -lopenblas
	LIBFLAGS_TEST := -undefined dynamic_lookup -lpthread -lm -lomp -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/openblas/lib -lopenblas
	
	EXTRA_LIBFLAGS := -L$(LIBRARY_PATH) -l$(EXTLIBNAME)
else
	ifeq ($(UNAME_OS),Linux)
		ifdef SOURSE_BLAS_IS_INTEL
$(info compiling for Intel libs on Linux)
			CFLAGS_APP := -I$(SRCDIR) -I$(EXTSRC) -I$(EXTSRC2)
			CFLAGS_PYLIB := -DPYBIND $(PYBINDINCL) -I$(EXTSRC) -I$(EXTSRC2) -fPIC
			
			EXTRA_CFLAGS := -O3 -qopenmp -DUTEST -Dintelmkl -Wall -DMKL_ILP64 -I"${MKLROOT}/include" -std=c++17  -pthread
						
			LIBFLAGS_PY := -shared -qopenmp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
			LIBFLAGS_TEST := -qopenmp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
			
			EXTRA_LIBFLAGS := -L$(LIBRARY_PATH) -l$(EXTLIBNAME)

		else # using openBLAS
			CFLAGS_APP := -I$(SRCDIR) -I$(EXTSRC) -I$(EXTSRC2)
			CFLAGS_PYLIB := -DPYBIND $(PYBINDINCL) -I$(EXTSRC) -I$(EXTSRC2) -fPIC
			EXTRA_CFLAGS := -Wall -Wextra -pedantic -D UTEST -03 -fopenmp -I/opt/openBLAS/include -std=c++17  -pthread
			LIBFLAGS_PY := -shared -lpthread -lm -lgomp /opt/openBLAS/lib/libopenblas.a
			LIBFLAGS_TEST := -lpthread -lm -lgomp /opt/openBLAS/lib/libopenblas.a

			EXTRA_LIBFLAGS := -L$(LIBRARY_PATH) -l$(EXTLIBNAME)
		endif
	endif
endif

ifdef perf_num
	EXTRA_CFLAGS +=
endif

TARGET_PYLIB := $(BDIR_PYLIB)/$(APPNAME)$(PYSUFF)
TARGET_APP := $(BDIR_APP)/$(APPNAME)_app.exe
TARGET_TEST := $(BDIR_TST)/$(APPNAME)_test.exe

CXX: $(TARGET_PYLIB) $(TARGET_TEST) $(TARGET_APP)

#----------------------------------------------------------------
SOURCES_ALL := $(wildcard $(SRCDIR)/*.cpp)

SOURCES = $(filter-out $(SRCDIR)/evoped_app.cpp, $(wildcard $(SRCDIR)/*.cpp))

OBJECT_PYL := $(patsubst $(SRCDIR)/%,$(BDIR_PYLIB)/%,$(SOURCES:.cpp=.o))
OBJECT_PYL += $(BDIR_PYLIB)/sparse_matrix.o

OBJECT_APP := $(filter-out $(BDIR_APP)/wrapper.o, $(patsubst $(SRCDIR)/%,$(BDIR_APP)/%,$(SOURCES_ALL:.cpp=.o)))
OBJECT_APP += $(BDIR_APP)/sparse_matrix.o

SOURCES_TST := $(wildcard $(TESTDIR)/$(SRCDIR)/*.cpp)
OBJECT_TST := $(patsubst $(TSRCDIR)/%,$(BDIR_TST)/%,$(SOURCES_TST:.cpp=.o))
OBJECT_TST += $(patsubst $(SRCDIR)/%,$(BDIR_PYLIB)/%,$(SOURCES:.cpp=.o))
OBJECT_TST += $(BDIR_APP)/sparse_matrix.o

#----------------------------------------------------------------
# build lib for python
$(TARGET_PYLIB): $(OBJECT_PYL)
	@echo " Linking python lib ..."
	$(CXX) $^ $(LIBFLAGS_PY) -o $(TARGET_PYLIB) $(EXTRA_LIBFLAGS)

# build catch2 test
$(TARGET_TEST): $(OBJECT_TST)
	@echo " Linking tests ..."
	$(CXX) $^ $(LIBFLAGS_TEST) $(EXTRA_LIBFLAGS) -o $(TARGET_TEST)

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

# --------------------------------------------------------------
# Compiling external c++ file for APP
$(BDIR_APP)/%.o: $(EXTSRC)/sparse_matrix.cpp
	@echo " Compilling c++ app ..."
	@mkdir -p $(BDIR_APP)
	$(CXX) -c $(CFLAGS_APP) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling external c++ file for python lib
$(BDIR_PYLIB)/%.o: $(EXTSRC)/sparse_matrix.cpp
	@echo " Compilling pylib ..."
	@mkdir -p $(BDIR_PYLIB)
	$(CXX) -c $(CFLAGS_PYLIB) $(EXTRA_CFLAGS) -o $@ $< 

#----------------------------------------------------------------

clean:
	@echo " Cleaning..."
	@rm -fr $(BDIR_PYLIB) $(BDIR_TST) $(BDIR_APP) $(TARGET_PYLIB) $(TARGET_TEST) 2>/dev/null || true
#	@rm amatrix_* py_out_solution_model_* cpp_solution_model_* 2>/dev/null || true

test:
	@echo " Testing ..."
#	python3 $(TESTDIR)/$(SRCDIR)/*.py
ifeq ($(UNAME_OS),Darwin)
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
