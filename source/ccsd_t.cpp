#include "ccsd_t.hpp"

void CCSD_T::initialize()
{
    t3d.resize( nocc, nocc, nocc, nvirt, nvirt, nvirt );
    t3c.resize( nocc, nocc, nocc, nvirt, nvirt, nvirt );
    Dijkabc.resize( nocc, nocc, nocc, nvirt, nvirt, nvirt );
}

// P(i/jk) * P(a/bc)
// (ijk - jik - kji) * (abc - bac - cba)
// (ijkabc - ijkbac - ijkcba - jikabc + jikbac + jikcba - kjiabc + kjibac + kjicba)
void CCSD_T::build_disconnected_triples()
{
    // i, j, k -- occupied
    // a, b, c -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int k = 0; k < nocc; ++k )
            {
                for ( int a = nocc; a < size_; ++a )
                {
                    for ( int b = nocc; b < size_; ++b )
                    {
                        for ( int c = nocc; c < size_; ++c )
                        {
                            double res = 0;
                            // ijkabc
                            res += t1(i, a-nocc) * utilities.ASTwoElectronMOIntegrals(j, k, b, c);
                            // ijkbac
                            res -= t1(i, b-nocc) * utilities.ASTwoElectronMOIntegrals(j, k, a, c);
                            // ijkcba
                            res -= t1(i, c-nocc) * utilities.ASTwoElectronMOIntegrals(j, k, b, a);
                            // jikabc
                            res -= t1(j, a-nocc) * utilities.ASTwoElectronMOIntegrals(i, k, b, c);
                            // jikbac
                            res += t1(j, b-nocc) * utilities.ASTwoElectronMOIntegrals(i, k, a, c);
                            // jikcba
                            res += t1(j, c-nocc) * utilities.ASTwoElectronMOIntegrals(i, k, b, a);
                            // kjiabc
                            res -= t1(k, a-nocc) * utilities.ASTwoElectronMOIntegrals(j, i, b, c);
                            // kjibac
                            res += t1(k, b-nocc) * utilities.ASTwoElectronMOIntegrals(j, i, a, c);
                            // kjicba
                            res += t1(k, c-nocc) * utilities.ASTwoElectronMOIntegrals(j, i, b, a);

                            t3d(i,j,k,a-nocc, b-nocc,c-nocc) = res;
                        }
                    }
                }
            }
        }
    }
}

// P(i/jk) * P(a/bc)
// (ijk - jik - kji) * (abc - bac - cba)
// (ijkabc - ijkbac - ijkcba - jikabc + jikbac + jikcba - kjiabc + kjibac + kjicba)
void CCSD_T::build_connected_triples()
{
    // i, j, k -- occupied
    // a, b, c -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int k = 0; k < nocc; ++k )
            {
                for ( int a = nocc; a < size_; ++a )
                {
                    for ( int b = nocc; b < size_; ++b )
                    {
                        for ( int c = nocc; c < size_; ++c )
                        {
                            double sum1_ijkabc = 0.0;
                            double sum1_ijkbac = 0.0;
                            double sum1_ijkcba = 0.0;
                            double sum1_jikabc = 0.0;
                            double sum1_jikbac = 0.0;
                            double sum1_jikcba = 0.0;
                            double sum1_kjiabc = 0.0;
                            double sum1_kjibac = 0.0;
                            double sum1_kjicba = 0.0;

                            for ( int e = nocc; e < size_; ++e )
                            {
                                sum1_ijkabc += t2(j, k, a-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,i,b,c);
                                sum1_ijkbac -= t2(j, k, b-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,i,a,c);
                                sum1_ijkcba -= t2(j, k, c-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,i,b,a);
                                sum1_jikabc -= t2(i, k, a-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,j,b,c);
                                sum1_jikbac += t2(i, k, b-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,j,a,c);
                                sum1_jikcba += t2(i, k, c-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,j,b,a);
                                sum1_kjiabc -= t2(j, i, a-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,k,b,c);
                                sum1_kjibac += t2(j, i, b-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,k,a,c);
                                sum1_kjicba += t2(j, i, c-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(e,k,b,a);
                            }

                            double sum2_ijkabc = 0.0;
                            double sum2_ijkbac = 0.0;
                            double sum2_ijkcba = 0.0;
                            double sum2_jikabc = 0.0;
                            double sum2_jikbac = 0.0;
                            double sum2_jikcba = 0.0;
                            double sum2_kjiabc = 0.0;
                            double sum2_kjibac = 0.0;
                            double sum2_kjicba = 0.0;
                            for ( int m = 0; m < nocc; ++m )
                            {
                                sum2_ijkabc += t2(i, m, b-nocc, c-nocc) * utilities.ASTwoElectronMOIntegrals(m,a,j,k);
                                sum2_ijkbac -= t2(i, m, a-nocc, c-nocc) * utilities.ASTwoElectronMOIntegrals(m,b,j,k);
                                sum2_ijkcba -= t2(i, m, b-nocc, a-nocc) * utilities.ASTwoElectronMOIntegrals(m,c,j,k);
                                sum2_jikabc -= t2(j, m, b-nocc, c-nocc) * utilities.ASTwoElectronMOIntegrals(m,a,i,k);
                                sum2_jikbac += t2(j, m, a-nocc, c-nocc) * utilities.ASTwoElectronMOIntegrals(m,b,i,k);
                                sum2_jikcba += t2(j, m, b-nocc, a-nocc) * utilities.ASTwoElectronMOIntegrals(m,c,i,k);
                                sum2_kjiabc -= t2(k, m, b-nocc, c-nocc) * utilities.ASTwoElectronMOIntegrals(m,a,j,i);
                                sum2_kjibac += t2(k, m, a-nocc, c-nocc) * utilities.ASTwoElectronMOIntegrals(m,b,j,i);
                                sum2_kjicba += t2(k, m, b-nocc, a-nocc) * utilities.ASTwoElectronMOIntegrals(m,c,j,i);
                            }

                            t3c(i, j, k, a-nocc, b-nocc, c-nocc) = sum1_ijkabc - sum2_ijkabc + \
                                                                   sum1_ijkbac - sum2_ijkbac + \
                                                                   sum1_ijkcba - sum2_ijkcba + \
                                                                   sum1_jikabc - sum2_jikabc + \
                                                                   sum1_jikbac - sum2_jikbac + \
                                                                   sum1_jikcba - sum2_jikcba + \
                                                                   sum1_kjiabc - sum2_kjiabc + \
                                                                   sum1_kjibac - sum2_kjibac + \
                                                                   sum1_kjicba - sum2_kjicba;
                        }
                    }
                }
            }
        }
    }
}

void CCSD_T::build_Dijkabc()
{
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int k = 0; k < nocc; ++k )
            {
                for ( int a = nocc; a < size_; ++a )
                {
                    for ( int b = nocc; b < size_; ++b )
                    {
                        for ( int c = nocc; c < size_; ++c )
                        {
                            //std::cout << "Fock(c,c): " << utilities.SOFockMatrix(c, c) << std::endl;
                            Dijkabc(i,j,k,a-nocc,b-nocc,c-nocc) = utilities.SOFockMatrix(i, i) + \
                                                                  utilities.SOFockMatrix(j, j) + \
                                                                  utilities.SOFockMatrix(k, k) - \
                                                                  utilities.SOFockMatrix(a, a) - \
                                                                  utilities.SOFockMatrix(b, b) - \
                                                                  utilities.SOFockMatrix(c, c);
                            //std::cout << "i: " << i << "; j: " << j << "; a: " << a << "; b: " << b << "; c: " << c << " " << Dijkabc(i,j,k,a,b,c) << std::endl;
                        }
                    }
                }
            }
        }
    }
}

double CCSD_T::compute_perturbation()
{
    double res = 0.0;
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int k = 0; k < nocc; ++k )
            {
                for ( int a = nocc; a < size_; ++a )
                {
                    for ( int b = nocc; b < size_; ++b )
                    {
                        for ( int c = nocc; c < size_; ++c )
                        {
                            res += t3c(i,j,k,a-nocc,b-nocc,c-nocc) * \
                                   (t3c(i,j,k,a-nocc,b-nocc,c-nocc) + t3d(i,j,k,a-nocc,b-nocc,c-nocc)) / \
                                    Dijkabc(i,j,k,a-nocc,b-nocc,c-nocc);
                        }
                    }
                }
            }
        }
    }

    return res / 36.0;
}
