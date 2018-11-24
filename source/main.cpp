#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "molecule.hpp"

#include "ccsd.hpp"

int main()
{
    std::cout << std::fixed << std::setprecision(12);
    std::clock_t start = std::clock();

	Basis basis;
    //basis.read("./basis/h2o_cc_pvdz.gamess-us.dat");
    //basis.read("./basis/h2o_sto3g_gamess-us.dat");
    basis.read("./basis/sto3g-second-period-gamess-us.dat");
    basis.show("short");

    Molecule molecule;
    molecule.setBasis( &basis );

    molecule.readGeometryFile("geometry/h2o_crawford.dat");
    //molecule.readGeometryFile("geometry/ch4_crawford.dat");
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
    double SCF_energy = molecule.SCF();
    std::cout << "SCF Energy (total): " << SCF_energy << std::endl;

    // MP2
    molecule.fillTwoElectronMOIntegrals();
    double MP2_correction = molecule.computeMP2_correction();
    std::cout << "(main) MP2_correction: " << MP2_correction << std::endl;

    // CCSD
    int size_ = 2 * molecule.size(); // количество спинорбиталей
    int nocc = molecule.get_charge(); // количество занятых спинорбиталей
    int nvirt = size_ - nocc; // количество свободных (виртуальных) спинорбиталей
    CCSD ccsd( size_, nocc, nvirt );

    ccsd.initialize();
    ccsd.preparation( molecule );

    ccsd.run_diis();

    //double CCSD_correction = ccsd.run();
    //std::cout << "Total CCSD correction: " << CCSD_correction << std::endl;

    std::cout << "Total time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

    return 0;
}


