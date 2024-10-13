APPNAME := plinkio

SRCDIR := src
BDIR_CLIB := release_lib
BDIR_APP := release_app

# determine OS and source of BLAS lib
UNAME_OS := $(shell uname -s)
SOURSE_BLAS_IS_INTEL := 0

# determine compiler
ifeq ($(UNAME_OS),Darwin)
	CXX := clang++
else
#	ifeq ( $(SOURSE_BLAS_IS_INTEL), 1)
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
	CFLAGS_APP := -I$(SRCDIR)
	CFLAGS_CLIB := -I$(SRCDIR)

	EXTRA_CFLAGS := -O3 -Xpreprocessor -std=c++17 -Wall -Wextra -pedantic -Wno-c99-extensions

	LIBFLAGS := -undefined dynamic_lookup -shared -lm -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
	LIBFLAGS_TEST := -undefined dynamic_lookup -lm -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
	EXTRA_LIBFLAGS :=  -L$(BDIR_CLIB) -l$(APPNAME)
else
	ifeq ($(UNAME_OS),Linux)
		ifdef SOURSE_BLAS_IS_INTEL
			CFLAGS_APP := -O3 -I$(SRCDIR) -DMKL_ILP64 -I${MKLROOT}/include -std=c++17
			CFLAGS_CLIB := -O3 -I$(SRCDIR) -DMKL_ILP64 -I${MKLROOT}/include -fPIC -std=c++17
			
			EXTRA_CFLAGS := -Wall
			EXTRA_LIBFLAGS := -L$(BDIR_CLIB) -l$(APPNAME)
			EXTRA_LIBFLAGS_2 :=
			
			LIBFLAGS := -shared -lm
			LIBFLAGS_TEST := -lm
		else # using openBLAS
			CFLAGS_APP := -I$(SRCDIR) -I$(EXTSRC)
			EXTRA_CFLAGS := -Wall -Wextra -pedantic -03 -std=c++17
			LIBFLAGS_TEST := -lm
		endif
	endif
endif

ifdef perf_num
	EXTRA_CFLAGS += -D PERF_TEST
endif

TARGET_CLIB := $(BDIR_CLIB)/lib$(APPNAME)$(VERSUFF)
TARGET_APP := $(BDIR_APP)/$(APPNAME)_app.exe

CXX: $(TARGET_CLIB) $(TARGET_APP)

SOURCES_ALL := $(wildcard $(SRCDIR)/*.cpp)
SOURCES = $(filter-out $(SRCDIR)/plinkio_app.cpp, $(wildcard $(SRCDIR)/*.cpp))
OBJECT_LIB := $(patsubst $(SRCDIR)/%,$(BDIR_CLIB)/%,$(SOURCES:.cpp=.o))
OBJECT_APP := $(patsubst $(SRCDIR)/%,$(BDIR_APP)/%,$(SOURCES_ALL:.cpp=.o))

#----------------------------------------------------------------

# build shared c++ lib
$(TARGET_CLIB): $(OBJECT_LIB)
	@echo " Linking shared lib ..."
	$(CXX) $^ $(LIBFLAGS) -o $(TARGET_CLIB)

# build c++ app
$(TARGET_APP): $(OBJECT_APP)
	@echo " Linking c++ app ..."
	$(CXX) $^ $(LIBFLAGS_TEST) $(EXTRA_LIBFLAGS) -o $(TARGET_APP)

#----------------------------------------------------------------

# Compiling c++ lib
$(BDIR_CLIB)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling clib ..."
	@mkdir -p $(BDIR_CLIB)
	$(CXX) -c $(CFLAGS_CLIB) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling c++ app
$(BDIR_APP)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling c++ app ..."
	@mkdir -p $(BDIR_APP)
	$(CXX) -c $(CFLAGS_APP) $(EXTRA_CFLAGS) -o $@ $< 

#----------------------------------------------------------------

clean:
	@echo " Cleaning..."
	@rm -fr $(BDIR_CLIB) $(BDIR_APP) 2>/dev/null || true
	@rm amatrix_* py_out_solution_model_* cpp_solution_model_* 2>/dev/null || true

rclean:
	@echo " Cleaning..."
	@rm amatrix_* py_out_solution_model_* cpp_solution_model_* 2>/dev/null || true

run:
	@echo " Running ..."
	@export  LD_LIBRARY_PATH=$(PWDADDR)/$(BDIR_CLIB):"$(LIBRARY_PATH)"; $(TARGET_APP) test_data/YY_plink

cleanlogs:
	@echo " Cleaning logs ..."
	@rm *.log

.PHONY: clean
