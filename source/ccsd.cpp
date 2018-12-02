#include "ccsd.hpp"

void CCSD::initialize()
{
    Dia.resize( nocc, nvirt );
    Dijab.resize( nocc, nocc, nvirt, nvirt );

    t1.resize( nocc, nvirt );
    t2.resize( nocc, nocc, nvirt, nvirt );

    t1_updated.resize( nocc, nvirt );
    t2_updated.resize( nocc, nocc, nvirt, nvirt );

    tau.resize( nocc, nocc, nvirt, nvirt );
    tilda_tau.resize( nocc, nocc, nvirt, nvirt );

    Fae.resize( nvirt, nvirt );
    Fme.resize( nocc, nvirt );
    Fmi.resize( nocc, nocc );

    Wmnij.resize( nocc, nocc, nocc, nocc );
    Wabef.resize( nvirt, nvirt, nvirt, nvirt );
    Wmbej.resize( nocc, nvirt, nvirt, nocc );
}

void CCSD::fill_Dia()
{
    // a -- unoccupied, i -- occupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int a = nocc; a < size_; ++a )
        {
            Dia(i, a - nocc) = utilities.SOFockMatrix(i, i) - utilities.SOFockMatrix(a, a);
        }
    }
}

void CCSD::fill_Dijab()
{
    // i, j -- occupied
    // a, b -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = nocc; a < size_; ++a )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    Dijab(i, j, a - nocc, b - nocc) = utilities.SOFockMatrix(i, i) + utilities.SOFockMatrix(j, j) - \
                                                      utilities.SOFockMatrix(a, a) - utilities.SOFockMatrix(b, b);
                }
            }
        }
    }
}

void CCSD::fill_initial_t1()
{
    // i - occupied
    // a - unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int a = nocc; a < size_; ++a )
        {
            t1(i, a - nocc) = 0.0;
        }
    }
}

void CCSD::fill_initial_t2( Eigen::VectorXd const & HF_OrbitalEnergies )
{
    // i, j -- occupied
    // a, b -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = nocc; a < size_; ++a )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    // HF_OrbitalEnergies(i/2) -- потому что дважды занятые орбитали перешли в две спин-орбитали
                    t2(i, j, a - nocc, b - nocc) = utilities.ASTwoElectronMOIntegrals(i, j, a, b) / \
                            (HF_OrbitalEnergies(i/2) + HF_OrbitalEnergies(j/2) - HF_OrbitalEnergies(a/2) - HF_OrbitalEnergies(b/2));
                }
            }
        }
    }
}

double CCSD::test_MP2_Energy()
{
    double res = 0.0;

    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = nocc; a < size_; ++a )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    res += utilities.ASTwoElectronMOIntegrals(i, j, a, b) * t2(i, j, a - nocc, b - nocc);
                    //std::cout << "twoelectronMOIN: " << ASTwoElectronMOIntegrals(i, j, a, b) << std::endl;
                    //std::cout << "t2: " << t2(i, j, a, b) << std::endl;
                }
            }
        }
    }

    return res / 4.0;
}

void CCSD::fill_tau()
{
    // i, j -- occupied
    // a, b -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = 0; a < nvirt; ++a )
            {
                for( int b = 0; b < nvirt; ++b )
                {
                    tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a) * t1(j, b) -  t1(i, b) * t1(j, a);
                }
            }
        }
    }
}

void CCSD::fill_tilda_tau()
{
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = 0; a < nvirt; ++a )
            {
                for ( int b = 0; b < nvirt; ++b )
                {
                    tilda_tau(i, j, a, b) = t2(i, j, a, b) + 0.5 * (t1(i, a) * t1(j, b) - t1(i, b) * t1(j, a));
                }
            }
        }
    }
}

void CCSD::fill_Fae()
{
    // a, e -- unoccupied
    for ( int a = nocc; a < size_; ++a )
    {
        for ( int e = nocc; e < size_; ++e )
        {
            double sum1 = 0.0;
            // m -- occupied
            for ( int m = 0; m < nocc; ++m )
            {
                sum1 += utilities.SOFockMatrix(m, e) * t1(m, a - nocc);
            }
            sum1 *= 0.5;

            double sum2 = 0.0;
            // m -- occupied
            // f -- unoccupied
            for ( int m = 0; m < nocc; ++m )
            {
                for ( int f = nocc; f < size_; ++f )
                {
                    sum2 += t1(m, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, a, f, e);
                }
            }

            double sum3 = 0.0;
            // m,n -- occupied
            // f -- unoccupied
            for ( int m = 0; m < nocc; ++m )
            {
                for ( int n = 0; n < nocc; ++n )
                {
                    for ( int f = nocc; f < size_; ++f )
                    {
                        sum3 += tilda_tau(m, n, a-nocc, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, e, f);
                    }
                }
            }
            sum3 *= 0.5;

            Fae(a-nocc, e-nocc) = (1-delta(a, e)) * utilities.SOFockMatrix(a, e) - sum1 + sum2 - sum3;
        }
    }
}

void CCSD::fill_Fme()
{
    // m -- occupied
    // e -- unoccupied
    for ( int m = 0; m < nocc; ++m )
    {
        for ( int e = nocc; e < size_; ++e )
        {
            double sum1 = 0.0;
            // n -- occupied
            // f -- unoccupied
            for ( int n = 0; n < nocc; ++n )
            {
                for ( int f = nocc; f < size_; ++f )
                {
                    sum1 += t1(n, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, e, f);
                }
            }

            Fme(m, e-nocc) = utilities.SOFockMatrix(m, e) + sum1;
        }
    }
}

void CCSD::fill_Fmi()
{
    // i, m -- occupied
    for ( int m = 0; m < nocc; ++m )
    {
        for ( int i = 0; i < nocc; ++i )
        {
            double sum1 = 0.0;
            // e -- unoccupied
            for ( int e = nocc; e < size_; ++e )
            {
                sum1 += t1(i, e-nocc) * utilities.SOFockMatrix(m, e);
            }
            sum1 *= 0.5;

            double sum2 = 0.0;
            // e -- unoccupied
            // n -- occupied
            for ( int e = nocc; e < size_; ++e )
            {
                for ( int n = 0; n < nocc; ++n )
                {
                    sum2 += t1(n, e-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, i, e);
                }
            }

            double sum3 = 0.0;
            // n -- occupied
            // e, f -- unoccupied
            for ( int e = nocc; e < size_; ++e )
            {
                for ( int f = nocc; f < size_; ++f )
                {
                    for ( int n = 0; n < nocc; ++n )
                    {
                        sum3 += tilda_tau(i, n, e-nocc, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, e, f);
                    }
                }
            }
            sum3 *= 0.5;

            Fmi(m, i) = (1-delta(m, i)) * utilities.SOFockMatrix(m, i) + sum1 + sum2 + sum3;
        }
    }
}

void CCSD::fill_Wmnij()
{
    // i, j, m, n -- occupied
    for ( int m = 0; m < nocc; ++m )
    {
        for ( int n = 0; n < nocc; ++n )
        {
            for ( int i = 0; i < nocc; ++i )
            {
                for ( int j = 0; j < nocc; ++j )
                {
                    double sum1 = 0.0;
                    // e -- unoccupied
                    for ( int e = nocc; e < size_; ++e )
                    {
                        sum1 += t1(j, e-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, i, e);
                    }

                    double sum2 = 0.0;
                    // e -- unoccupied
                    for ( int e = nocc; e < size_; ++e )
                    {
                        sum2 += t1(i, e-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, j, e);
                    }

                    double sum3 = 0.0;
                    // e, f -- unoccupied
                    for ( int e = nocc; e < size_; ++e )
                    {
                        for ( int f = nocc; f < size_; ++f )
                        {
                            sum3 += tau(i, j, e-nocc, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, e, f);
                        }
                    }
                    sum3 *= 0.25;

                    Wmnij(m, n, i, j) = utilities.ASTwoElectronMOIntegrals(m, n, i, j) + sum1 - sum2 + sum3;
                }
            }
        }
    }
}

void CCSD::fill_Wabef()
{
    // a, b, e, f -- unoccupied
    for ( int a = nocc; a < size_; ++a )
    {
        for ( int b = nocc; b < size_; ++b )
        {
            for ( int e = nocc; e < size_; ++e )
            {
                for ( int f = nocc; f < size_; ++f )
                {
                    double sum1 = 0.0;
                    // m -- occupied
                    for ( int m = 0; m < nocc; ++m )
                    {
                        sum1 += t1(m, b-nocc) * utilities.ASTwoElectronMOIntegrals(a, m, e, f);
                    }

                    double sum2 = 0.0;
                    // m -- occupied
                    for ( int m = 0; m < nocc; ++m )
                    {
                        sum2 += t1(m, a-nocc) * utilities.ASTwoElectronMOIntegrals(b, m, e, f);
                    }

                    double sum3 = 0.0;
                    // m, n -- occupied
                    for ( int m = 0; m < nocc; ++m )
                    {
                        for ( int n = 0; n < nocc; ++n )
                        {
                            sum3 += tau(m, n, a-nocc, b-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, e, f);
                        }
                    }
                    sum3 *= 0.25;

                    Wabef(a-nocc, b-nocc, e-nocc, f-nocc) = utilities.ASTwoElectronMOIntegrals(a, b, e, f) - sum1 + sum2 + sum3;
                }
            }
        }
    }
}

void CCSD::fill_Wmbej()
{
    // m, j -- occupied
    // b, e -- unccupied
    for ( int m = 0; m < nocc; ++m )
    {
        for ( int b = nocc; b < size_; ++b )
        {
            for ( int e = nocc; e < size_; ++e )
            {
                for ( int j = 0; j < nocc; ++j )
                {
                    double sum1 = 0.0;
                    // f -- unoccupied
                    for ( int f = nocc; f < size_; ++f )
                    {
                        sum1 += t1(j, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, b, e, f);
                    }

                    double sum2 = 0.0;
                    // n -- occupied
                    for ( int n = 0; n < nocc; ++n )
                    {
                        sum2 += t1(n, b-nocc) * utilities.ASTwoElectronMOIntegrals(m, n, e, j);
                    }

                    double sum3 = 0.0;
                    // n -- occupied
                    // f -- unoccupied
                    for ( int n = 0; n < nocc; ++n )
                    {
                        for ( int f = nocc; f < size_; ++f )
                        {
                            sum3 += (0.5 * t2(j, n, f-nocc, b-nocc) + t1(j, f-nocc) * t1(n, b-nocc)) * \
                                    utilities.ASTwoElectronMOIntegrals(m, n, e, f);
                        }
                    }

                    Wmbej(m, b-nocc, e-nocc, j) = utilities.ASTwoElectronMOIntegrals(m, b, e, j) + sum1 - sum2 - sum3;
                }
            }
        }
    }
}

void CCSD::update_t1()
{
    // i -- occupied
    // a -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int a = nocc; a < size_; ++a )
        {
            double sum1 = 0.0;
            // e -- unoccupied
            for ( int e = nocc; e < size_; ++e )
            {
                sum1 += t1(i, e-nocc) * Fae(a-nocc, e-nocc);
            }

            double sum2 = 0.0;
            // m -- occupied
            for ( int m = 0; m < nocc; ++m )
            {
                sum2 += t1(m, a-nocc) * Fmi(m, i);
            }

            double sum3 = 0.0;
            // m -- occupied, e -- unoccupied
            for ( int m = 0; m < nocc; ++m )
            {
                for ( int e = nocc; e < size_; ++e )
                {
                    sum3 += t2(i, m, a-nocc, e-nocc) * Fme(m, e-nocc);
                }
            }

            double sum4 = 0.0;
            // n -- occupied, f -- unoccupied
            for ( int n = 0; n < nocc; ++n )
            {
                for ( int f = nocc; f < size_; ++f )
                {
                    sum4 += t1(n, f-nocc) * utilities.ASTwoElectronMOIntegrals(n, a, i, f);
                }
            }

            double sum5 = 0.0;
            // m -- occupied; e, f -- unoccupied
            for ( int m = 0; m < nocc; ++m )
            {
                for ( int e = nocc; e < size_; ++e )
                {
                    for ( int f = nocc; f < size_; ++f )
                    {
                        sum5 += t2(i, m, e-nocc, f-nocc) * utilities.ASTwoElectronMOIntegrals(m, a, e, f);
                    }
                }
            }
            sum5 *= 0.5;

            double sum6 = 0.0;
            // m, n -- occupied; e -- unoccupied
            for ( int m = 0; m < nocc; ++m )
            {
                for ( int n = 0; n < nocc; ++n )
                {
                    for ( int e = nocc; e < size_; ++e )
                    {
                        sum6 += t2(m, n, a-nocc, e-nocc) * utilities.ASTwoElectronMOIntegrals(n, m, e, i);
                    }
                }
            }
            sum6 *= 0.5;

            t1_updated(i, a-nocc) = utilities.SOFockMatrix(i, a) + sum1 - sum2 + sum3 - sum4 - sum5 - sum6;
            t1_updated(i, a-nocc) /= Dia(i, a-nocc);
        }
    }
}

void CCSD::update_t2()
{
    // i, j -- occupied; a, b -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = nocc; a < size_; ++a )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    double sum1ab = 0.0;
                    // меняем местами a и b (действие оператора P_(a, b)
                    double sum1ba = 0.0;

                    // e -- unoccupied
                    for ( int e = nocc; e < size_; ++e )
                    {
                        double res = 0.0;
                        // m -- occupied
                        for ( int m = 0; m < nocc; ++m )
                        {
                            res += t1(m, b-nocc) * Fme(m, e-nocc);
                        }
                        sum1ab += t2(i, j, a-nocc, e-nocc) * (Fae(b-nocc, e-nocc) - 0.5 * res);

                        res = 0.0;
                        // m -- occupied
                        for ( int m = 0; m < nocc; ++m )
                        {
                            res += t1(m, a-nocc) * Fme(m, e-nocc);
                        }
                        sum1ba += t2(i, j, b-nocc, e-nocc) * (Fae(a-nocc, e-nocc) - 0.5 * res);
                    }

                    double sum2ij = 0.0;
                    // меняем местами i и j (действие оператора P_(i, j)
                    double sum2ji = 0.0;

                    // m -- occupied
                    for ( int m = 0; m < nocc; ++m )
                    {
                        double res = 0.0;
                        // e -- unoccupied
                        for ( int e = nocc; e < size_; ++e )
                        {
                            res += t1(j, e-nocc) * Fme(m, e-nocc);
                        }

                        sum2ij += t2(i, m, a-nocc, b-nocc) * (Fmi(m, j) + 0.5 * res);

                        res = 0.0;
                        for ( int e = nocc; e < size_; ++e )
                        {
                            res += t1(i, e-nocc) * Fme(m, e-nocc);
                        }

                        sum2ji += t2(j, m, a-nocc, b-nocc) * (Fmi(m, i) + 0.5 * res);
                    }

                    double sum3 = 0.0;
                    // m, n -- occupied
                    for ( int m = 0; m < nocc; ++m )
                    {
                        for ( int n = 0; n < nocc; ++n )
                        {
                            sum3 += tau(m, n, a-nocc, b-nocc) * Wmnij(m, n, i, j);
                        }
                    }
                    sum3 *= 0.5;

                    double sum4 = 0.0;
                    // e, f -- unoccupied
                    for ( int e = nocc; e < size_; ++e )
                    {
                        for ( int f = nocc; f < size_; ++f )
                        {
                            sum4 += tau(i, j, e-nocc, f-nocc) * Wabef(a-nocc, b-nocc, e-nocc, f-nocc);
                        }
                    }
                    sum4 *= 0.5;

                    // дальше у нас идут 4 вида сумм из-за того, что стоит произведение операторов
                    // P_(i, j) * P_(a, b).
                    // Сначала идет самый простой вариант, когда никакие индексы не переставлены
                    double sum5ijab = 0.0;
                    // m -- occupied, e -- unoccupied
                    // Переставляем местами b и a
                    double sum5ijba = 0.0;
                    // Переставляем местами i и j
                    double sum5jiab = 0.0;
                    // Одновременно переставляем местами (i, j) и (a, b)
                    double sum5jiba = 0.0;

                    for ( int m = 0; m < nocc; ++m )
                    {
                        for ( int e = nocc; e < size_; ++e )
                        {
                            sum5ijab += (t2(i, m, a-nocc, e-nocc) * Wmbej(m, b-nocc, e-nocc, j) - \
                                         t1(i, e-nocc) * t1(m, a-nocc) * utilities.ASTwoElectronMOIntegrals(m, b, e, j));
                            sum5ijba += (t2(i, m, b-nocc, e-nocc) * Wmbej(m, a-nocc, e-nocc, j) - \
                                         t1(i, e-nocc) * t1(m, b-nocc) * utilities.ASTwoElectronMOIntegrals(m, a, e, j));
                            sum5jiab += (t2(j, m, a-nocc, e-nocc) * Wmbej(m, b-nocc, e-nocc, i) - \
                                         t1(j, e-nocc) * t1(m, a-nocc) * utilities.ASTwoElectronMOIntegrals(m, b, e, i));
                            sum5jiba += (t2(j, m, b-nocc, e-nocc) * Wmbej(m, a-nocc, e-nocc, i) - \
                                         t1(j, e-nocc) * t1(m, b-nocc) * utilities.ASTwoElectronMOIntegrals(m, a, e, i));
                        }
                    }

                    // 2 суммы, оператор P_(i, j)
                    double sum6ij = 0.0;
                    // переставляем местами (i, j)
                    double sum6ji = 0.0;

                    // e -- unoccupied
                    for ( int e = nocc; e < size_; ++e )
                    {
                        sum6ij += (t1(i, e-nocc) * utilities.ASTwoElectronMOIntegrals(a, b, e, j));
                        sum6ji += (t1(j, e-nocc) * utilities.ASTwoElectronMOIntegrals(a, b, e, i));
                    }

                    // 2 суммы, оператор P_(a, b)
                    double sum7ab = 0.0;
                    // переставляем местами (a, b)
                    double sum7ba = 0.0;

                    // m -- occupied
                    for ( int m = 0; m < nocc; ++m )
                    {
                        sum7ab += (t1(m, a-nocc) * utilities.ASTwoElectronMOIntegrals(m, b, i, j));
                        sum7ba += (t1(m, b-nocc) * utilities.ASTwoElectronMOIntegrals(m, a, i, j));
                    }

                    t2_updated(i, j, a-nocc, b-nocc) = utilities.ASTwoElectronMOIntegrals(i, j, a, b) + sum1ab - sum1ba - sum2ij + sum2ji + \
                        sum3 + sum4 + sum5ijab - sum5ijba - sum5jiab + sum5jiba + sum6ij - sum6ji - sum7ab + sum7ba;
                    t2_updated(i, j, a-nocc, b-nocc) /= Dijab(i, j, a-nocc, b-nocc);
                }
            }
        }
    }
}

double CCSD::computeCCSD_correction()
{
    double sum1 = 0.0;
    // i -- occupied, a -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int a = nocc; a < size_; ++a )
        {
            sum1 += utilities.SOFockMatrix(i, a) * t1_updated(i, a-nocc);
        }
    }
    //std::cout << "(CCSD correction) sum1: " << sum1 << std::endl;

    double sum2 = 0.0;
    // i, j -- occupied; a, b -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = nocc; a < size_; ++a )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    sum2 += utilities.ASTwoElectronMOIntegrals(i, j, a, b) * t2_updated(i, j, a-nocc, b-nocc);
                }
            }
        }
    }
    sum2 *= 0.25;
    //std::cout << "(CCSD correction) sum2: " << sum2 << std::endl;

    double sum3 = 0.0;
    // i, j -- occupied; a, b -- unoccupied
    for ( int i = 0; i < nocc; ++i )
    {
        for ( int j = 0; j < nocc; ++j )
        {
            for ( int a = nocc; a < size_; ++a )
            {
                for ( int b = nocc; b < size_; ++b )
                {
                    sum3 += utilities.ASTwoElectronMOIntegrals(i, j, a, b) * t1_updated(i, a-nocc) * t1_updated(j, b-nocc);
                }
            }
        }
    }
    sum3 *= 0.5;
    //std::cout << "(CCSD correction) sum3: " << sum3 << std::endl;

    return sum1 + sum2 + sum3;
}

void CCSD::preparation( Molecule const & molecule, MP2 const& mp2 )
{
    utilities.fillAS_MO_TwoElectronIntegrals( mp2.get_two_electron_MO_integrals() );
    utilities.fillSOHcore( molecule.get_C(), molecule.get_Hcore() );
    utilities.fillSOFock();

    fill_Dia();
    //std::cout << "(main) Dia: " << std::endl << ccsd.get_Dia() << std::endl;

    fill_Dijab();
    //std::cout << "(main) Diab: " << std::endl << ccsd.get_Dijab() << std::endl;

    // initial guess for t1 and t2 cluster amplitudes
    fill_initial_t1();
    fill_initial_t2( molecule.get_HF_OrbitalEnergies() );
    double test_MP2_correction = test_MP2_Energy();
    std::cout << "(preparation) testing MP2 correcion: " << test_MP2_correction << std::endl;

    t1_updated = t1;
    t2_updated = t2;
}

void CCSD::iterate()
{
    fill_tau();
    fill_tilda_tau();

    fill_Fae();
    fill_Fme();
    fill_Fmi();

    fill_Wmnij();
    fill_Wabef();
    fill_Wmbej();

    update_t1();
    update_t2();
}

double CCSD::run()
{
    double CCSDcorr_E_old = 0.0;
    double CCSDcorr_E = 0.0;

    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Starting CCSD iterations" << std::endl;

    auto ccsd_start = std::chrono::high_resolution_clock::now();

    for ( int CCSD_iter = 0; CCSD_iter < maxiter; ++CCSD_iter, CCSDcorr_E_old = CCSDcorr_E )
    {
        iterate();

        CCSDcorr_E = computeCCSD_correction();

        std::cout << "CCSD iteration: " << CCSD_iter << "; CCSD correlation: " << CCSDcorr_E <<
                     "; DE: " << CCSDcorr_E - CCSDcorr_E_old << std::endl;

        if ( std::abs(CCSDcorr_E - CCSDcorr_E_old) < E_conv )
            break;

        // обновляем кластерные операторы
        t1 = t1_updated;
        t2 = t2_updated;
    }

    auto ccsd_end = std::chrono::high_resolution_clock::now();
    std::cout << "CCSD iterations took " << std::chrono::duration_cast<std::chrono::milliseconds>( ccsd_end - ccsd_start ).count() / 1000.0
              << " s." << std::endl;

    return CCSDcorr_E;
}

double CCSD::find_max( Eigen::MatrixXd const & B )
{
    double max_ = std::abs(B(0, 0));

    for ( int i = 0; i < B.rows() - 1; ++i )
    {
        for ( int j = 0; j < B.cols() - 1; ++j )
        {
            if ( std::abs(B(i, j)) > max_ )
                max_ = std::abs(B(i, j));
        }
    }

    return max_;
}

double CCSD::run_diis()
{
    auto ccsd_diis_start = std::chrono::high_resolution_clock::now();

    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Starting CCSD iterations (DIIS acceleration)" << std::endl;

    std::vector<Eigen::MatrixXd> diis_vals_t1;
    std::vector<Eigen::Tensor<double, 4>> diis_vals_t2;
    diis_vals_t1.push_back( t1 );
    diis_vals_t2.push_back( t2 );

    double CCSDcorr_E = computeCCSD_correction();
    double CCSDcorr_E_old = CCSDcorr_E;

    std::cout << "CCSD iteration: 0; CCSD correlation: " << CCSDcorr_E <<
                 "; DE: " << CCSDcorr_E - CCSDcorr_E_old << "; MP2: " << std::endl;

    size_t diis_size = 0;
    Eigen::MatrixXd oldt1;
    Eigen::Tensor<double, 4> oldt2;

    std::vector<Eigen::VectorXd> diis_errors;

    for ( int CCSD_iter = 1; CCSD_iter < maxiter; ++CCSD_iter, CCSDcorr_E_old = CCSDcorr_E )
    {
        oldt1 = t1;
        oldt2 = t2;

        //std::cout << "oldt1: \n" << oldt1 << std::endl;

        iterate();
        t1 = t1_updated;
        t2 = t2_updated;

        CCSDcorr_E = computeCCSD_correction();

        std::cout << "CCSD iteration: " << CCSD_iter << "; CCSD correlation: " << CCSDcorr_E <<
                     "; DE: " << CCSDcorr_E - CCSDcorr_E_old << "; DIIS = " << diis_size << std::endl;

        if ( std::abs(CCSDcorr_E - CCSDcorr_E_old) < E_conv_diis )
            break;

        diis_vals_t1.push_back( t1 );
        diis_vals_t2.push_back( t2 );

        // build error vector
        Eigen::MatrixXd diff1 = t1 - oldt1;
        // пакуем эту разницу по строкам в Map
        Eigen::Map<Eigen::VectorXd, Eigen::RowMajor> error_t1( diff1.data(), diff1.size() );
        //std::cout << "error_t1.size: " << error_t1.size() << std::endl;

        Eigen::Tensor<double, 4> diff2 = t2 - oldt2;
        // пакуем эту разницу по строкам в Map
        Eigen::Map<Eigen::VectorXd, Eigen::RowMajor> error_t2( diff2.data(), diff2.size() );
        //std::cout << "error_t2.size: " << error_t2.size() << std::endl;

        Eigen::VectorXd diis_error = Eigen::VectorXd::Zero( error_t1.size() + error_t2.size() );
        for ( int i = 0; i < error_t1.size(); ++i )
            diis_error(i) = error_t1(i);
        for ( int i = error_t1.size(); i < diis_error.size(); ++i )
            diis_error(i) = error_t2(i - error_t1.size());

        diis_errors.push_back(diis_error);

        // обрезаем diis_vals_t1, diis_vals_t2
        if ( static_cast<int>(diis_vals_t1.size()) > max_diis )
        {
            //std::cout << "Max diis size is reached. Truncating diis vector." << std::endl;

            diis_vals_t1.erase( diis_vals_t1.begin() );
            diis_vals_t2.erase( diis_vals_t2.begin() );
            diis_errors.erase( diis_errors.begin() );
        }

        diis_size = diis_vals_t1.size() - 1;
        //std::cout << "diis_size: " << diis_size << std::endl;

        // инициализируем матрицу DIIS минус единицами
        Eigen::MatrixXd B = Eigen::MatrixXd::Constant( diis_size + 1, diis_size + 1, -1.0 );
        B( diis_size, diis_size ) = 0.0;

        for ( size_t i = 0; i < diis_errors.size(); ++i )
        {
            for ( size_t j = 0; j < diis_errors.size(); ++j )
            {
                double res = 0.0;

                for ( int k = 0; k < diis_errors[i].size(); ++k )
                    res += diis_errors[i](k) * diis_errors[j](k);

                B(i, j) = res;
            }
        }

        double max_ = find_max( B );

        for ( int i = 0; i < B.rows() - 1; ++i )
            for ( int j = 0; j < B.cols() - 1; ++j )
                B(i, j) /= max_;

        //std::cout << "B: \n" << B << std::endl;

        Eigen::VectorXd resid = Eigen::VectorXd::Zero( diis_size + 1 );
        resid( diis_size ) = -1.0;

        Eigen::VectorXd ci = B.colPivHouseholderQr().solve( resid );

        //std::cout << "ci: \n" << ci << std::endl;

        t1.setZero();
        t2.setZero();
        for ( size_t i = 0;  i < diis_size; ++i )
        {
            t1 += ci(i) * diis_vals_t1[i + 1];
            t2 += ci(i) * diis_vals_t2[i + 1];
        }

        //std::cout << "t1:\n" << t1 << std::endl;
    }

    auto ccsd_diis_end = std::chrono::high_resolution_clock::now();

    std::cout << "------------------------------------------" << std::endl;
    std::cout << "CCSD (DIIS) iterations took " << std::chrono::duration_cast<std::chrono::milliseconds>( ccsd_diis_end - ccsd_diis_start ).count() / 1000.0
              << " s." << std::endl;

    return CCSDcorr_E;
}
