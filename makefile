C = g++
CFLAGS = -g --std=c++11
EIGEN = -I /usr/local/include/eigen3


INPUT = input.dat
OUTPUT = basis.log
DIR = result

cc-pvDZ = ./tests/basis/cc-pvdz.gamess-us.dat
cc-pvTZ = ./tests/basis/cc-pvtz.gamess-us.dat

BASIS = ./source/Nikita/basis.cpp
GEOM = ./source/Nikita/geom.cpp

basis_EXE = ./basis
geom_EXE = ./geom
FILES = *.log *.dat

all: clean test

test: basis test1 test2 MKDIR

	@$(basis_EXE) < $(INPUT) > $(OUTPUT)
	@mv $(FILES) $(basis_EXE) $(DIR)

MKDIR:

	@mkdir $(DIR)

test1: 

	@echo $(cc-pvDZ) >> $(INPUT)

test2: 

	@echo $(cc-pvTZ) >> $(INPUT)

.PHONY: basis

basis: $(BASIS)

	@$(C) $(CFLAGS) $(BASIS) -o $(basis_EXE)

geom: $(GEOM)

	@$(C) $(CFLAGS) $(EIGEN) $(GEOM) -o $(geom_EXE)

.PHONY: clean

clean:

	@rm -rf $(DIR) $(FILES) $(basis_EXE) $(geom_EXE)
