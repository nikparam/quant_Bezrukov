.PHONY: clean all

CCXX = g++

CXXFLAGS = -std=c++11 -g -O2 -lm -Wall -Wextra
LINK_EIGEN = -I /usr/local/include/eigen3/

BUILDDIR := ./build/

all: main

main: $(addprefix $(BUILDDIR), $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
	$(CCXX) $(CXXFLAGS) $(LINK_EIGEN) $^ -o $@

$(BUILDDIR)%.o: %.cpp
	@echo ">> (g++) COMPILING $<...";
	@$(CCXX) $(CXXFLAGS) $(LINK_EIGEN) -c -MD $< -o $@

include $(wildcard *.d)

clean:
	@rm $(BUILDDIR)*.o
	@rm $(BUILDDIR)*.d
