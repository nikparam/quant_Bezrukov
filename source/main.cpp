#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "molecule.hpp"

int main()
{
    std::cout << std::fixed << std::setprecision(12);
    std::clock_t start = std::clock();

	Basis basis;
    // basis.read("./basis/h2o_cc_pvdz.gamess-us.dat");
    basis.read("./basis/h2o_sto3g_gamess-us.dat");
    basis.show("full");

    Molecule molecule;
    molecule.setBasis( &basis );

    molecule.readGeometryFile("geometry/h2o_crawford.dat");
    molecule.showAtoms();

    molecule.setCharge();
    molecule.setOutput("./out.txt");

    molecule.fillOverlapMatrix();
    molecule.fillKineticEnergyMatrix();
    molecule.fillNuclearAttractionMatrix();
    molecule.fillElectronRepulsionMatrix();
    molecule.showElectronAttractionMatrix();

    // создаем матрицу P
    molecule.makeInitialGuess();

    molecule.SCF_initialize();
    molecule.SCF();

    // MP2
    molecule.fillTwoElectronMOIntegrals();
    double MP2_correction = molecule.computeMP2_correction();
    std::cout << "(main) MP2_correction: " << MP2_correction << std::endl;

    // CCSD
    molecule.fillAS_MO_TwoElectronIntegrals();
    molecule.fillSOHcoreMatrix();
    molecule.fillSOFockMatrix();

    std::cout << "(main) SO F: " << std::endl << molecule.SOFockMatrix << std::endl;

    // initial guess for tia and tijab cluster amplitudes
    molecule.fill_initial_tia();
    molecule.fill_initial_tijab();
    double test_MP2_correction = molecule.testCCSD_MP2_Energy();
    std::cout << "(main) testing MP2 correcion: " << test_MP2_correction << std::endl;

    molecule.fill_tau_ijab();
    molecule.fill_tau_tilda_ijab();
    molecule.fill_D1();
    molecule.fill_D2();

    molecule.fill_F();
    std::cout << "(main) F: " << std::endl << molecule.intermF << std::endl;

    molecule.fill_W();

    molecule.update_t1();
    molecule.update_t2();

    double CCSD_correction = molecule.computeCCSD_correction();
    std::cout << "(main) CCSD_correction: " << CCSD_correction << std::endl;

    std::cout << "Total time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

	return 0;
}


