#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"
#include "Basis.hpp"

int main()
{
	Basis basis;
	basis.read("./basis/cc-pvdz.gamess-us.dat");

    basis.show("short");

	return 0;
}


