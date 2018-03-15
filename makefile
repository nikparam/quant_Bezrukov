C = g++
CFLAGS = -g --std=c++11

INPUT = input.dat
OUTPUT = basis.log
DIR = result

cc-pvDZ = ./tests/basis/cc-pvdz.gamess-us.dat
cc-pvTZ = ./tests/basis/cc-pvtz.gamess-us.dat

CODE = ./source/Nikita/basis.cpp
EXECUTABLE = ./basis
FILES = *.log *.dat

all: clean test

test: basis test1 test2 MKDIR

	@$(EXECUTABLE) < $(INPUT) > $(OUTPUT)
	@mv $(FILES) $(EXECUTABLE) $(DIR)

MKDIR:

	@mkdir $(DIR)

test1: 

	@echo $(cc-pvDZ) >> $(INPUT)

test2: 

	@echo $(cc-pvTZ) >> $(INPUT)

.PHONY: basis

basis: $(CODE)

	@$(C) $(CFLAGS) $(CODE) -o $(EXECUTABLE)

.PHONY: clean

clean:

	@rm -rf $(DIR)
