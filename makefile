C = g++
CFLAGS = -g --std=c++11

basis: ./source/Nikita/basis.cpp

	$(C) $(CFLAGS) $^ -o basis

test1: ./tests/basis/cc-pvdz.gamess-us.dat

	echo '$^' > basis.log
	echo '$^' | ./basis >> basis.log

test2: ./tests/basis/cc-pvtz.gamess-us.dat

	echo '$^' > basis.log
	echo '$^' | ./basis >> basis.log

clean:
	rm -rf basis *.log
