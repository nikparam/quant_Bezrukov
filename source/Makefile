.PHONY: clean all

CCXX = g++

CXXFLAGS = -std=c++11 -g -O2 -lm -Wall -Wextra

BUILDDIR := ./build/

all: main

main: $(addprefix $(BUILDDIR), $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
	$(CCXX) $(CXXFLAGS) $^ -o $@

$(BUILDDIR)%.o: %.cpp
	@echo ">> (g++) COMPILING $<...";
	@$(CCXX) $(CXXFLAGS) -c -MD $< -o $@

include $(wildcard *.d)

clean:
	@rm $(BUILDDIR)*.o
	@rm $(BUILDDIR)*.d
