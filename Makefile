# Template for Athena++ Makefile
# The 'configure.py' python script uses this template to create the actual Makefile

# Files for conditional compilation

COORDINATES_FILE = spherical_polar.cpp
EOS_FILE = adiabatic_mhd.cpp
PROBLEM_FILE = disk.cpp
RECONSTRUCT_FILE = plm.cpp
RSOLVER_FILE = hlld.cpp
RSOLVER_DIR = mhd/
HYDRO_INT_FILE = vl2.cpp

# General compiler specifications

CXX := g++
CPPFLAGS := 
CXXFLAGS := -O3
LDFLAGS := 
LDLIBS := 

# Preliminary definitions

EXE_DIR := bin/
EXECUTABLE := $(EXE_DIR)athena
SRC_FILES := $(wildcard src/*.cpp) \
	     $(wildcard src/bvals/*.cpp) \
	     $(wildcard src/coordinates/*.cpp) \
	     src/eos/$(EOS_FILE) \
	     $(wildcard src/field/*.cpp) \
	     $(wildcard src/hydro/*.cpp) \
	     $(wildcard src/hydro/srcterms/*.cpp) \
	     src/hydro/rsolvers/$(RSOLVER_DIR)$(RSOLVER_FILE) \
	     $(wildcard src/mesh/*.cpp) \
	     $(wildcard src/outputs/*.cpp) \
	     $(wildcard src/reconstruct/*.cpp) \
	     $(wildcard src/task_list/*.cpp) \
	     $(wildcard src/utils/*.cpp) \
	     src/pgen/$(PROBLEM_FILE) \
	     src/pgen/default_pgen.cpp
OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
SRC_DIR := $(dir $(SRC_FILES) $(PROB_FILES))
VPATH := $(SRC_DIR)

# Generally useful targets

.PHONY : all dirs clean

all : dirs $(EXECUTABLE)

dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Link objects into executable

$(EXECUTABLE) : $(OBJ_FILES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $(OBJ_FILES) $(LDFLAGS) $(LDLIBS)

# Create objects from source files

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Cleanup

clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)
