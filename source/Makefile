.PHONY: clean all

F := gfortran
CCXX := g++

CXXFLAGS := -std=c++11 -g -O2 -lm -Wall -Wextra
LINK_EIGEN := -isystem /usr/local/include/eigen3/
FORFLAGS := -ffixed-line-length-none -ffree-line-length-none
LINKFLAGS := -lgfortran

BUILDDIR := ./build/

all: main

main: $(addprefix $(BUILDDIR), $(patsubst %.cpp, %.o, $(wildcard *.cpp))) \
	  $(addprefix $(BUILDDIR), $(patsubst %.f90, %.oo,$(wildcard *.f90)))
	$(CCXX) $(CXXFLAGS) $(LINK_EIGEN) $^ -o $@ $(LINKFLAGS)

$(BUILDDIR)%.o: %.cpp
	@echo ">> (g++) COMPILING $<...";
	@$(CCXX) $(CXXFLAGS) $(LINK_EIGEN) -c -MD $< -o $@

$(BUILDDIR)%.oo: %.f90
	@echo ">> (gfortran) COMPILING STUPID FORTRAN $<...";
	@$(F) $(FORFLAGS) -c $< -o $@

include $(wildcard *.d)

clean:
	@rm $(BUILDDIR)*.o
	@rm $(BUILDDIR)*.d
