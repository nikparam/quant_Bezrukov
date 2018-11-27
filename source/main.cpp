#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "molecule.hpp"

#include "ccsd_utilities.hpp"
#include "ccsd.hpp"
#include "ccsd_t.hpp"
#include "cis.hpp"

int main()
{
    std::cout << std::fixed << std::setprecision(12);
    std::clock_t start = std::clock();

    Basis basis;
    //basis.DEBUG = 1;

    //basis.read("./basis/h2o_cc_pvdz.gamess-us.dat");
    //basis.read("./basis/h2o_sto3g_gamess-us.dat");
    basis.read("./basis/sto3g-second-period-gamess-us.dat");
    //basis.read("./basis/basis_6-31G.txt");
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
    //molecule.showElectronAttractionMatrix();

    /*
    Eigen::Tensor<double, 4> eri = molecule.get_two_electron_integrals();
    for ( int i = 0; i < eri.dimension(0); ++i )
    {
        for ( int j = 0; j < eri.dimension(1); ++j )
        {
            std::cout << eri(i, j, 0, 0) << " ";
        }
        std::cout << std::endl;
    }
    */

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

    // Вспомогательный класс, хранящий:
    // SOHcore, SOFock, антисимметризованные двуэлектронные интегралы на молекулярных орбиталях
    CCSD_Utilities ccsd_utilities( size_, nocc, nvirt );
    // эта последовательность выполняется в ccsd.prepation()
    //ccsd_utilities.fillAS_MO_TwoElectronIntegrals( molecule.get_two_electron_MO_integrals() );
    //ccsd_utilities.fillSOHcore( molecule.get_C(), molecule.get_Hcore() );
    //ccsd_utilities.fillSOFock();

    // CIS
    /*
    CIS cis( size_, nocc, nvirt, ccsd_utilities );
    cis.initialize();
    cis.fill_cis_matrix();

    Eigen::MatrixXd & cis_matrix = cis.get_cis_matrix();
    for ( int i = 0; i < cis_matrix.rows(); ++i )
    {
        for ( int j = 0; j < cis_matrix.cols(); ++j )
        {
            std::cout << cis_matrix(i,j) << " ";
        }
        std::cout << std::endl;
    }
    */

    /*
    Eigen::VectorXd eigs = cis.diagonalize_cis_matrix();
    std::cout << "CIS matrix eigenvalues: " << std::endl;
    for ( int k = 0; k < eigs.size(); ++k )
        std::cout << eigs(k) << std::endl;
    */

    CCSD ccsd( size_, nocc, nvirt, ccsd_utilities );
    ccsd.initialize(); // memory allocation
    ccsd.preparation( molecule );

    //double CCSD_correction = ccsd.run_diis();
    //std::cout << "Total CCSD correction (with DIIS acceleration): " << CCSD_correction << std::endl;

    double CCSD_correction = ccsd.run();
    std::cout << "Total CCSD correction: " << CCSD_correction << std::endl;

    double totalCCSD_energy = SCF_energy + CCSD_correction;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Total CCSD energy: " << totalCCSD_energy << std::endl;

    // CCSD(T)
    CCSD_T ccsd_t( size_, nocc, nvirt, ccsd_utilities );
    ccsd_t.initialize();
    ccsd_t.set_t1( ccsd.get_t1_updated() );
    ccsd_t.set_t2( ccsd.get_t2_updated() );

    ccsd_t.build_disconnected_triples();
    ccsd_t.build_connected_triples();
    ccsd_t.build_Dijkabc();

    double CCSD_T_correction = ccsd_t.compute_perturbation();
    std::cout << "CCSD(T) correction: " << CCSD_T_correction << std::endl;

    double totalCCSD_T_energy = SCF_energy + CCSD_correction + CCSD_T_correction;
    std::cout << "Total CCSD(T) energy: " << totalCCSD_T_energy << std::endl;

    std::cout << "Total time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

    //Eigen::Tensor<double, 6> & t3d = ccsd_t.get_t3d();
    //for ( int i = 0; i < t3d.dimension(0); ++i )
    //{
    //    for ( int j = 0; j < t3d.dimension(1); ++j )
    //    {
    //        std::cout << t3d(i, j, 1, 2, 0, 3) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    return 0;
}


