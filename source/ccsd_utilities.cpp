#include "ccsd_utilities.hpp"

void CCSD_Utilities::initialize_utilities()
{
    // количество спинорбиталей вдвое больше количества дважды занятых орбиталей
    ASTwoElectronMOIntegrals.resize( size_, size_, size_, size_ );
    SOHcoreMatrix.resize( size_, size_ );
    SOFockMatrix.resize( size_, size_ );
}

void CCSD_Utilities::fillAS_MO_TwoElectronIntegrals( Eigen::Tensor<double, 4> const & twoElectronMOIntegrals )
// Antisymmetrized integrals over spin-orbitals
{
    double int1, int2;

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            for ( int r = 0; r < size_; ++r )
            {
                for ( int s = 0; s < size_; ++s )
                {
                    int1 = twoElectronMOIntegrals(p/2, r/2, q/2, s/2) * (r%2 == p%2) * (s%2 == q%2);
                    int2 = twoElectronMOIntegrals(p/2, s/2, q/2, r/2) * (p%2 == s%2) * (q%2 == r%2);
                    ASTwoElectronMOIntegrals(p, q, r, s) = int1 - int2;
                }
            }
        }
    }
}

void CCSD_Utilities::fillSOHcore( Eigen::MatrixXd const & matrixC, Eigen::MatrixXd const & matrixHcore )
{
    Eigen::MatrixXd tmp = matrixC.transpose() * matrixHcore * matrixC;

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            SOHcoreMatrix(p, q) = tmp(p/2, q/2) * (p%2 == q%2);
        }
    }
}

void CCSD_Utilities::fillSOFock()
// fill spin-orbital Fock matrix
{
    double res = 0;

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            res = SOHcoreMatrix(p, q);
            for ( int m = 0; m < nocc; ++m )
                res += ASTwoElectronMOIntegrals(p, m, q, m);

            SOFockMatrix(p, q) = res;
        }
    }
}



