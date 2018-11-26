#include "cis.hpp"

void CIS::initialize()
{
    cis_matrix.resize( nvirt * nocc, nvirt * nocc );
}

void CIS::fill_cis_matrix()
{
    int ii = 0;
    int jj = 0;
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int a = nocc; a < size_; ++a )
        {
            for ( int j = 0; j < nocc; ++j )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    cis_matrix(ii,jj) += utilities.SOFockMatrix(a,b) * delta(i,j) - \
                                         utilities.SOFockMatrix(i,j) * delta(a,b) + \
                                         utilities.ASTwoElectronMOIntegrals(a,j,i,b);
                    //std::cout << "ii: " << ii << "; jj: " << jj << "; i: " << i << "; a: " << a << "; j: " << j << "; b: " << b << std::endl;

                    ++jj;
                    if ( jj == cis_matrix.rows() )
                        jj = 0;
                }
            }
            ++ii;
        }
    }
}

Eigen::VectorXd CIS::diagonalize_cis_matrix()
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es( cis_matrix );
    return es.eigenvalues();
}
