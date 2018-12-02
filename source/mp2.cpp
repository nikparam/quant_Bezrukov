#include "mp2.hpp"

void MP2::fillTwoElectronMOIntegrals_eff()
{
    Eigen::Tensor<double, 4> X1, X2, X3; // промежуточное хранение
    X1.resize( size_, size_, size_, size_ );
    X2.resize( size_, size_, size_, size_ );
    X3.resize( size_, size_, size_, size_ );
    twoElectronMOIntegrals.resize( size_, size_, size_, size_);

    // самый внутренний цикл
    for ( int mu = 0; mu < size_; ++mu )
    {
        for ( int nu = 0; nu < size_; ++nu )
        {
            for ( int lambda = 0; lambda < size_; ++lambda )
            {
                for ( int s = 0; s < size_; ++s )
                {
                    double res = 0;
                    for ( int sigma = 0; sigma < size_; ++sigma )
                        res += matrixC(sigma, s) * electronRepulsionTensor(mu, nu, lambda, sigma);
                    X1(mu, nu, lambda, s) = res;
                }
            }
        }
    }

    for ( int mu = 0; mu < size_; ++mu )
    {
        for ( int nu = 0; nu < size_; ++nu )
        {
            for ( int s = 0; s < size_; ++s )
            {
                for ( int r = 0; r < size_; ++r )
                {
                    double res = 0.0;
                    for ( int lambda = 0; lambda < size_; ++lambda )
                        res += matrixC(lambda, r) * X1(mu, nu, lambda, s);
                    X2(mu, nu, r, s) = res;
                }
            }
        }
    }

    for ( int mu = 0; mu < size_; ++mu )
    {
        for ( int q = 0; q < size_; ++q )
        {
            for ( int s = 0; s < size_; ++s )
            {
                for ( int r = 0; r < size_; ++r )
                {
                    double res = 0.0;
                    for ( int nu = 0; nu < size_; ++nu )
                        res += matrixC(nu, q) * X2(mu, nu, r, s);
                    X3(mu, q, r, s) = res;
                }
            }
        }
    }

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            for ( int s = 0; s < size_; ++s )
            {
                for ( int r = 0; r < size_; ++r )
                {
                    double res = 0.0;
                    for ( int mu = 0; mu < size_; ++mu )
                        res += matrixC(mu, p) * X3(mu, q, r, s);

                    twoElectronMOIntegrals(p, q, r, s) = res;
                }
            }
        }
    }
}


void MP2::fillTwoElectronMOIntegrals()
{
    twoElectronMOIntegrals.resize( size_, size_, size_, size_);

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            for ( int r = 0; r < size_; ++r )
            {
                for ( int s = 0; s < size_; ++s )
                {
                    double res = 0.0;

                    for ( int mu = 0; mu < size_; ++mu )
                    {
                        for ( int nu = 0; nu < size_; ++nu )
                        {
                            for ( int lambda = 0; lambda < size_; ++lambda )
                            {
                                for ( int sigma = 0; sigma < size_; ++sigma )
                                {
                                    res += matrixC(mu, p) * matrixC(nu, q) * matrixC(lambda, r) * matrixC(sigma, s) \
                                            * electronRepulsionTensor(mu, nu, lambda, sigma);
                                }
                            }
                        }
                    }

                    twoElectronMOIntegrals(p, q, r, s) = res;
                }
            }
        }
    }
}

double MP2::computeMP2_correction()
{
    double num = 0.0;
    double MP2_Energy = 0.0;

    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    num = twoElectronMOIntegrals(i, a, j, b) * (2.0 * twoElectronMOIntegrals(i, a, j, b) - \
                                                                twoElectronMOIntegrals(i, b, j, a));
                    MP2_Energy += num / (HF_OrbitalEnergies(i) + HF_OrbitalEnergies(j) - HF_OrbitalEnergies(a) - \
                                         HF_OrbitalEnergies(b));
                }
            }
        }
    }

    return MP2_Energy;
}
