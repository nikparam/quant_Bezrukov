#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "molecule.h"

int main()
{
	Basis basis;
    // basis.read("./basis/cc-pvdz.gamess-us.dat");
    basis.read("./basis/h2_cc-pvdz.gamess-us.dat");

    // basis.show("short");

    Molecule molecule;
    molecule.setBasis( &basis );

    molecule.readGeometryFile("geometry/h2.dat");
    molecule.showAtoms();

    molecule.fillOverlapMatrix();
    molecule.saveOverlapMatrix();

	return 0;
}


