#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "molecule.hpp"

int main()
{
	Basis basis;
    // basis.read("./basis/h2_cc-pvtz.gamess-us.dat");
    basis.read("./basis/h2o_sto3g_gamess-us.dat");
    // basis.read("./basis/h2_cc-pvdz.gamess-us.dat");
    basis.show("full");

    Molecule molecule;
    molecule.setBasis( &basis );

    molecule.readGeometryFile("geometry/h2.dat");
    molecule.showAtoms();

    molecule.fillOverlapMatrix();
    molecule.saveOverlapMatrix();

    molecule.fillKineticEnergyMatrix();
    molecule.saveKineticEnergyMatrix();

    molecule.fillNuclearAttractionMatrix();
    molecule.saveNuclearAttractionMatrix();

    molecule.fillElectronRepulsionMatrix();
    molecule.showElectronAttractionMatrix();

	return 0;
}


