C = g++
CFLAGS = -g --std=c++11

cc-pvDZ = './tests/basis/cc-pvdz.gamess-us.dat'
cc-pvTZ = './tests/basis/cc-pvtz.gamess-us.dat'

CODE = ./source/Nikita/basis.cpp
EXECUTABLE = basis
FILES = *.log

test1: $(EXECUTABLE)

	echo $(cc-pvDZ) > basis.log
	echo $(cc-pvDZ) | $(EXECUTABLE) >> basis.log

test2: $(EXECUTABLE)

	echo $(cc-pvTZ) > basis.log
	echo $(cc-pvTZ) | $(EXECUTABLE) >> basis.log


basis: $(CODE)

	$(C) $(CFLAGS) $(CODE) -o $(EXECUTABLE)

clean:

	rm -rf $(EXECUTABLE) $(FILES)
