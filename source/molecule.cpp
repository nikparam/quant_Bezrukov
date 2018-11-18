#include "molecule.hpp"

Molecule::Molecule()
{
}

void Molecule::readGeometryFile( std::string const & filename )
{
   std::ifstream infile( filename );
   if ( !infile )
       throw std::invalid_argument("Can't open the geometry file.");
   else
       parseGeometryFile( infile );

   infile.close();
}

void Molecule::parseGeometryFile( std::ifstream & infile )
{
    char buf[MAXLINE];
    std::string current_line;

    // пропускаем первые два элемента
    infile.getline(buf, MAXLINE);
    infile.getline(buf, MAXLINE);

    std::stringstream ss;
    std::string name;
    double x, y, z;

    while( infile.getline(buf, MAXLINE) )
    {
        current_line = buf;

        // почему-то просто ss.str("") не очищает буфер stringstream
        // а вот с ss.clear() работает ?!
        ss.str("");
        ss.clear();

        ss << current_line;
        ss >> name >> x >> y >> z;

        atoms.emplace_back(name, x, y, z);
        atoms.end()[-1].attachElementFunctions( basis );
    }
}

double Molecule::calculateEijt( int i, int j, int t, double Qx, double a, double b )
/*
    Recursive definition of Hermite Gaussian coefficients.
        Returns double.
        a: orbital exponent on Gaussian 'a'
        b: orbital exponent on Gaussian 'b'
        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
        t: number of nodes in Hermite (depends on type of integral,
           e.g. always zero for overlap integrals)
        Qx: distance between origins of Gaussian 'a' and 'b'
 */
{
    double p = a + b;
    double q = a * b / p;

    if ( t < 0 || t > i + j )
        return 0;

    else if ( (i == j) && (j == t) && (t == 0) )
        return std::exp( -q * Qx * Qx ); // K_AB

    else if ( j == 0 )
        // decrement index i
        return (1.0/(2.0 * p)) * calculateEijt(i - 1, j, t - 1, Qx, a, b) - \
                (q * Qx / a) * calculateEijt(i - 1, j, t, Qx, a, b) + \
                (t + 1) * calculateEijt(i - 1, j, t + 1, Qx, a, b);
    else
        // decrement index j
        return (1.0/(2.0 * p)) * calculateEijt(i, j - 1, t - 1, Qx, a, b) + \
                (q * Qx / b) * calculateEijt(i, j - 1, t, Qx, a, b) + \
                (t + 1) * calculateEijt(i, j - 1, t + 1, Qx, a, b);
}

double Molecule::overlapPrimitive( Primitive * a, Primitive * b,
                                   QuantumNumbers * t1, QuantumNumbers * t2,
                                   std::vector<double> coords1,
                                   std::vector<double> coords2 )
/*
Evaluates overlap integral between two Gaussians
  Returns a double
  a:   primitive: coefficient and exponent
  b:   primitive: coefficient and exponent
  t1:  three unsigned ints -- orbital angular momentum (e.g. (1,0,0) for Gaussian 'a'
  t2:  three unsigned ints -- orbital angular momentum for Gaussian 'b'
  coords1: a vector containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
  coords2: a vector containing origin of Gaussian 'b'
*/
{
    double S1 = calculateEijt(t1->i, t2->i, 0, coords1[0] - coords2[0], a->get_exponent(), b->get_exponent());
    double S2 = calculateEijt(t1->j, t2->j, 0, coords1[1] - coords2[1], a->get_exponent(), b->get_exponent());
    double S3 = calculateEijt(t1->k, t2->k, 0, coords1[2] - coords2[2], a->get_exponent(), b->get_exponent());
    return S1 * S2 * S3 * std::pow(M_PI / (a->get_exponent() + b->get_exponent()), 1.5);
}

double Molecule::overlapCGO( ContractedGaussianOrbital * cgo1, ContractedGaussianOrbital * cgo2,
                             std::vector<double> coords1, std::vector<double> coords2 )
{
    double res = 0.0;
    for ( size_t i = 0; i < cgo1->getPrimitivesCount(); i++ )
    {
        // коэффициент перед первым примитивом
        double coeff1 = cgo1->getPrimitive(i).get_coefficient();

        // тройка квантовых чисел
        QuantumNumbers t1 = cgo1->getQuantumNumbers();

        for ( size_t j = 0; j < cgo2->getPrimitivesCount(); j++ )
        {
            // коэффициент перед вторым примитивом
            double coeff2 = cgo2->getPrimitive(j).get_coefficient();

            // тройка квантовых чисел
            QuantumNumbers t2 = cgo2->getQuantumNumbers();

            res += coeff1 * coeff2 * overlapPrimitive( &cgo1->getPrimitive(i), &cgo2->getPrimitive(j),
                                                      &t1, &t2, coords1, coords2 );
        }
    }

    return res;
}

void Molecule::fillOverlapMatrix( )
{
    int size_ = size();
    overlapMatrix.resize(size_, size_);

    size_t it1 = 0, it2 = 0;

    // первый итератор по атомам
    for ( size_t i = 0; i < atoms.size(); i++ )
    {
        // итератор по контрактированным орбиталям первого атома
        for ( size_t j = 0; j < atoms[i].get_element()->getCGOCount(); j++ )
        {
            // второй итератор по атомам
            for ( size_t k = 0; k < atoms.size(); k++ )
            {
                // итератор по контрактированным орбиталям второго атома
                for ( size_t l = 0; l < atoms[k].get_element()->getCGOCount(); l++ )
                {
                    overlapMatrix(it1, it2) = overlapCGO(atoms[i].get_element()->getCGO(j),
                                                         atoms[k].get_element()->getCGO(l),
                                                         atoms[i].getCoords(),
                                                         atoms[k].getCoords() );
                    ++it2;
                    if ( (int) it2 == size_ )
                    {
                        ++it1;
                        it2 = 0;
                    }

                }
            }
        }
    }

    outFile << "S matrix: " << std::endl << overlapMatrix << std::endl << std::endl;
}


double Molecule::kineticPrimitive(Primitive * a, Primitive * b,
                                  QuantumNumbers *t1, QuantumNumbers *t2,
                                  std::vector<double> coords1, std::vector<double> coords2)
// Evaluates kinetic energy integral betweeen two primitives
{
    QuantumNumbers t2i = t2->change_i(2);
    QuantumNumbers t2j = t2->change_j(2);
    QuantumNumbers t2k = t2->change_k(2);

    double term0 = b->get_exponent() * (2 * t2->getAngularMomentum() + 3) * \
            overlapPrimitive(a, b, t1, t2, coords1, coords2);

    double term1 = -2.0 * std::pow(b->get_exponent(), 2) * \
            ( overlapPrimitive(a, b, t1, &t2i, coords1, coords2) + \
              overlapPrimitive(a, b, t1, &t2j, coords1, coords2) + \
              overlapPrimitive(a, b, t1, &t2k, coords1, coords2) );

    t2i = t2->change_i(-2);
    t2j = t2->change_j(-2);
    t2k = t2->change_k(-2);
    double term2 = -0.5 * (t2->i * (t2->i - 1) * overlapPrimitive(a, b, t1, &t2i, coords1, coords2) + \
                           t2->j * (t2->j - 1) * overlapPrimitive(a, b, t1, &t2j, coords1, coords2) + \
                           t2->k * (t2->k - 1) * overlapPrimitive(a, b, t1, &t2k, coords1, coords2));

    return term0 + term1 + term2;
}

double Molecule::kineticCGO( ContractedGaussianOrbital *cgo1, ContractedGaussianOrbital *cgo2,
                              std::vector<double> coords1, std::vector<double> coords2 )
// evaluates kinetic energy between two contracted Gaussian functions
{
    double sum = 0.0;
    for ( size_t i = 0; i < cgo1->getPrimitivesCount(); ++i )
    {
        for ( size_t j = 0; j < cgo2->getPrimitivesCount(); ++j )
        {
            sum += cgo1->getPrimitive(i).get_coefficient() * cgo2->getPrimitive(j).get_coefficient() * \
                    kineticPrimitive( &cgo1->getPrimitive(i), &cgo2->getPrimitive(j),
                                      &cgo1->getQuantumNumbers(), &cgo2->getQuantumNumbers(),
                                      coords1, coords2 );
        }
    }

    return sum;
}

void Molecule::fillKineticEnergyMatrix( )
{
    //std::cout << "a1: " << a1->get_name() << "; " << a1->get_x() << " " << a1->get_y() << " " << a1->get_z() << std::endl;
    //std::cout << "a2: " << a2->get_name() << "; " << a2->get_x() << " " << a2->get_y() << " " << a2->get_z() << std::endl;

    int size_ = size();
    kineticEnergyMatrix.resize(size_, size_);

    size_t it1 = 0, it2 = 0;

    // первый итератор по атомам
    for ( size_t i = 0; i < atoms.size(); i++ )
    {
        // итератор по контрактированным орбиталям первого атома
        for ( size_t j = 0; j < atoms[i].get_element()->getCGOCount(); j++ )
        {
            // второй итератор по атомам
            for ( size_t k = 0; k < atoms.size(); k++ )
            {
                // итератор по контрактированным орбиталям второго атома
                for ( size_t l = 0; l < atoms[k].get_element()->getCGOCount(); l++ )
                {
                    kineticEnergyMatrix(it1, it2) = kineticCGO(atoms[i].get_element()->getCGO(j),
                                                               atoms[k].get_element()->getCGO(l),
                                                               atoms[i].getCoords(),
                                                               atoms[k].getCoords() );
                    ++it2;
                    if ( (int) it2 == size_ )
                    {
                        ++it1;
                        it2 = 0;
                    }

                }
            }
        }
    }

    outFile << "Kinetic energy matrix: " << std::endl << kineticEnergyMatrix << std::endl << std::endl;
}

void Molecule::showAtoms()
{
    for ( auto & atom : atoms )
    {
        std::cout << "Name: " << atom.get_name() << "; " <<
                                 atom.get_x() << " " <<
                                 atom.get_y() << " " <<
                                 atom.get_z() << std::endl;

        std::cout << "Element: " << atom.get_element()->getName() <<
                     "; CGOcount = " << atom.get_element()->getCGOCount() << std::endl;


    }
}

double Molecule::calculateHCintegral(int t, int u, int v,
                                    int n, double p, double PCx, double PCy, double PCz, double RPC)
// Returns the Coulomb auxiliary Hermite integral
// Arguments:
// 	t, u, v: order of Coulomb Hermite derivative in x, y, z (Helgaker, Taylor)
//  n: order of Boys function
//  PCx, PCy, PCz: Cartesian vector distance between gaussian composite center P and nuclear center C
//  RPC: distance between P and C
{
    double T = p * RPC * RPC;
    double val = 0.0;
    if ( (t == u) && (u == v) && (v == 0) )
        val += std::pow(- 2 * p, n) * MathUtils::BoysFunction(n, T);
    else if ( (t == u) && (u == 0) )
    {
        if ( v > 1 )
            val += (v - 1) * calculateHCintegral(t, u, v - 2, n + 1, p, PCx, PCy, PCz, RPC);
        val += PCz * calculateHCintegral(t, u, v - 1, n + 1, p, PCx, PCy, PCz, RPC);
    }
    else if ( t == 0 )
    {
        if ( u > 1 )
            val += (u - 1) * calculateHCintegral(t, u - 2, v, n + 1, p, PCx, PCy, PCz, RPC);
        val += PCy * calculateHCintegral(t, u - 1, v, n + 1, p, PCx, PCy, PCz, RPC);
    }
    else
    {
        if ( t > 1 )
            val += (t - 1) * calculateHCintegral(t - 2, u, v, n + 1, p, PCx, PCy, PCz, RPC);
        val += PCx * calculateHCintegral(t - 1, u, v, n + 1, p, PCx, PCy, PCz, RPC);
    }

    return val;
}

std::vector<double> Molecule::gaussianProductCenter( double a, std::vector<double> & A,
                                                     double b, std::vector<double> & B )
{
    std::vector<double> P;
    P.push_back( (a * A[0] + b * B[0]) / (a + b) );
    P.push_back( (a * A[1] + b * B[1]) / (a + b) );
    P.push_back( (a * A[2] + b * B[2]) / (a + b) );
    return P;
}

double Molecule::nuclearAttractionPrimitive(Primitive *a, QuantumNumbers *t1, std::vector<double> coords1,
                                            Primitive *b, QuantumNumbers *t2, std::vector<double> coords2,
                                            std::vector<double> nuclearCoords)
// Evaluates a nuclear matrix element between two gaussian primitives and one nuclear
{
    double p = a->get_exponent() + b->get_exponent();
    std::vector<double> P = gaussianProductCenter(a->get_exponent(), coords1, b->get_exponent(), coords2);
    double RPC = std::pow( std::pow(P[0] - nuclearCoords[0], 2) + \
                           std::pow(P[1] - nuclearCoords[1], 2) + \
                           std::pow(P[2] - nuclearCoords[2], 2), 0.5);

    //std::cout << "P: " << P[0] << " " << P[1] << " " << P[2] << std::endl;
    //std::cout << "RPC: " << RPC << std::endl;

    double val = 0.0;
    for ( int t = 0; t < (t1->i + t2->i + 1); ++t )
    {
        for ( int u = 0; u < (t1->j + t2->j + 1); ++u )
        {
            for ( int v = 0; v < (t1->k + t2->k + 1); ++v )
            {
                val += calculateEijt(t1->i, t2->i, t, coords1[0] - coords2[0], a->get_exponent(), b->get_exponent()) * \
                       calculateEijt(t1->j, t2->j, u, coords1[1] - coords2[1], a->get_exponent(), b->get_exponent()) * \
                       calculateEijt(t1->k, t2->k, v, coords1[2] - coords2[2], a->get_exponent(), b->get_exponent()) * \
                       calculateHCintegral(t, u, v, 0, p, P[0] - nuclearCoords[0], P[1] - nuclearCoords[1], P[2] - nuclearCoords[2], RPC);
            }
        }
    }

    val *= 2 * M_PI / p;
    return val;
}

double Molecule::nuclearAttractionCGO(ContractedGaussianOrbital *cgo1, ContractedGaussianOrbital *cgo2,
                                      std::vector<double> coords1, std::vector<double> coords2,
                                      std::vector<double> nuclearCoords )
// evaluates nuclear attraction between two contracted gaussians and a given nucleus
{
    double sum = 0.0;
    for ( size_t i = 0; i < cgo1->getPrimitivesCount(); ++i )
    {
        for ( size_t j = 0; j < cgo2->getPrimitivesCount(); ++j )
        {
            sum -= cgo1->getPrimitive(i).get_coefficient() * cgo2->getPrimitive(j).get_coefficient() * \
                    nuclearAttractionPrimitive( &cgo1->getPrimitive(i), &cgo1->getQuantumNumbers(), coords1,
                                                &cgo2->getPrimitive(j), &cgo2->getQuantumNumbers(), coords2,
                                                nuclearCoords );
        }
    }

    return sum;
}

double Molecule::nuclearAttractionCGO_total(ContractedGaussianOrbital *cgo1, ContractedGaussianOrbital *cgo2,
                                            std::vector<double> coords1, std::vector<double> coords2)
// evaluates nuclear attraction between two contracted gaussians and ALL nuclei
{
    double sum = 0.0;
    for ( size_t i = 0; i < atoms.size(); ++i )
    {
        //if ( !vectorsEqual(coords1, coords2, atoms[i].getCoords()) )
        //{
            sum += atoms[i].get_element()->getCharge() * \
                    nuclearAttractionCGO(cgo1, cgo2, coords1, coords2, atoms[i].getCoords() );
        //}
    }

    return sum;
}

void Molecule::fillNuclearAttractionMatrix( )
{
    int size_ = size();
    nuclearAttractionMatrix.resize(size_, size_);

    size_t it1 = 0, it2 = 0;

    // первый итератор по атомам
    for ( size_t i = 0; i < atoms.size(); i++ )
    {
        // итератор по контрактированным орбиталям первого атома
        for ( size_t j = 0; j < atoms[i].get_element()->getCGOCount(); j++ )
        {
            // второй итератор по атомам
            for ( size_t k = 0; k < atoms.size(); k++ )
            {
                // итератор по контрактированным орбиталям второго атома
                for ( size_t l = 0; l < atoms[k].get_element()->getCGOCount(); l++ )
                {
                    nuclearAttractionMatrix(it1, it2) = nuclearAttractionCGO_total( atoms[i].get_element()->getCGO(j),
                                                                                    atoms[k].get_element()->getCGO(l),
                                                                                    atoms[i].getCoords(),
                                                                                    atoms[k].getCoords() );
                    ++it2;
                    if ( (int) it2 == size_ )
                    {
                        ++it1;
                        it2 = 0;
                    }

                }
            }
        }
    }

    outFile << "Nuclear-Electron attraction matrix: " << std::endl << nuclearAttractionMatrix << std::endl << std::endl;
}

double Molecule::electronRepulsionPrimitive(Primitive *a, QuantumNumbers *ta, std::vector<double> A,
                                            Primitive *b, QuantumNumbers *tb, std::vector<double> B,
                                            Primitive *c, QuantumNumbers *tc, std::vector<double> C,
                                            Primitive *d, QuantumNumbers *td, std::vector<double> D)
// evaluates electron repulsion energy between primitives
{
    //std::ofstream out;
    //out.open("counter.txt", fstream::out | fstream::app);

    double p = a->get_exponent() + b->get_exponent();
    double q = c->get_exponent() + d->get_exponent();
    double alpha = p * q / (p + q);

    std::vector<double> P = gaussianProductCenter(a->get_exponent(), A, b->get_exponent(), B);
    std::vector<double> Q = gaussianProductCenter(c->get_exponent(), C, d->get_exponent(), D);

    double RPQ = std::pow( std::pow(P[0] - Q[0], 2) + std::pow(P[1] - Q[1], 2) + std::pow(P[2] - Q[2], 2), 0.5);

    double sum = 0.0;
    for ( int t = 0; t < (ta->i + tb->i + 1); ++t )
    {
        for ( int u = 0; u < (ta->j + tb->j + 1); ++u )
        {
            for ( int v = 0; v < (ta->k + tb->k + 1); ++v )
            {
                for ( int tau = 0; tau < (tc->i + td->i + 1); ++tau )
                {
                    for ( int nu = 0; nu < (tc->j + td->j + 1); ++nu )
                    {
                        for ( int phi = 0; phi < (tc->k + td->k + 1); ++phi )
                        {
                            sum += calculateEijt(ta->i, tb->i, t, A[0] - B[0], a->get_exponent(), b->get_exponent()) * \
                                   calculateEijt(ta->j, tb->j, u, A[1] - B[1], a->get_exponent(), b->get_exponent()) * \
                                   calculateEijt(ta->k, tb->k, v, A[2] - B[2], a->get_exponent(), b->get_exponent()) * \
                                   calculateEijt(tc->i, td->i, tau, C[0] - D[0], c->get_exponent(), d->get_exponent()) * \
                                   calculateEijt(tc->j, td->j, nu, C[1] - D[1], c->get_exponent(), d->get_exponent()) * \
                                   calculateEijt(tc->k, td->k, phi, C[2] - D[2], c->get_exponent(), d->get_exponent()) * \
                                   std::pow(-1.0, tau + nu + phi) * \
                                   calculateHCintegral(t + tau, u + nu, v + phi, 0, alpha, P[0] - Q[0], P[1] - Q[1], P[2] - Q[2], RPQ);
                        }
                    }
                }
            }
        }
    }

    sum *= 2.0 * std::pow(M_PI, 2.5) / (p * q * std::pow(p + q, 0.5));
    return sum;
}

double Molecule::electronRepulsionCGO(ContractedGaussianOrbital *a, std::vector<double> A,
                                      ContractedGaussianOrbital *b, std::vector<double> B,
                                      ContractedGaussianOrbital *c, std::vector<double> C,
                                      ContractedGaussianOrbital *d, std::vector<double> D )
// compute electron repulsion between CGOs
{
    double sum = 0.0;
    for ( size_t ja = 0; ja < a->getPrimitivesCount(); ++ja )
    {
        for ( size_t jb = 0; jb < b->getPrimitivesCount(); ++jb )
        {
            for ( size_t jc = 0; jc < c->getPrimitivesCount(); ++jc )
            {
                for ( size_t jd = 0; jd < d->getPrimitivesCount(); ++jd )
                {
                    sum += a->getPrimitive(ja).get_coefficient() * b->getPrimitive(jb).get_coefficient() * \
                           c->getPrimitive(jc).get_coefficient() * d->getPrimitive(jd).get_coefficient() * \
                           electronRepulsionPrimitive(&a->getPrimitive(ja), &a->getQuantumNumbers(), A,
                                                      &b->getPrimitive(jb), &b->getQuantumNumbers(), B,
                                                      &c->getPrimitive(jc), &c->getQuantumNumbers(), C,
                                                      &d->getPrimitive(jd), &d->getQuantumNumbers(), D );
                }
            }
        } 
   }

   return sum;
}

void Molecule::setOutput( std::string const & filename )
{
    std::string filename_ = filename;
    if ( filename_ == "" )
        filename_ = "out.txt";

    outFile.open(filename_);
    outFile << std::fixed << std::setprecision(8);
}

int Molecule::size( ) const
{
    int size = 0;
    for ( size_t i = 0; i < atoms.size(); i++ )
        size += atoms[i].get_element()->getCGOCount();

    return size;
}

void Molecule::fillElectronRepulsionMatrix( )
{
    electronRepulsionTensor.resize( size(), size(), size(), size() );

    int it1 = 0;
    // первый итератор по атомам
    for ( size_t a = 0; a < atoms.size(); ++a )
    {
        // итератор по контрактированным орбиталям первого атома
        for ( size_t i = 0; i < atoms[a].get_element()->getCGOCount(); ++i )
        {
            int it2 = 0;
            // второй итератор по атомам
            for ( size_t b = 0; b < atoms.size(); ++b )
            {
                // итератор по контрактированным орбиталям второго атома
                for ( size_t j = 0; j < atoms[b].get_element()->getCGOCount(); ++j )
                {
                    int it3 = 0;
                    // третий итератор по атомам
                    for ( size_t c = 0; c < atoms.size(); ++c )
                    {
                        // итератор по контрактированном орбиталям третьего атома
                        for ( size_t k = 0; k < atoms[c].get_element()->getCGOCount(); ++k )
                        {
                            int it4 = 0;
                            // четвертый итератор по атомам
                            for ( size_t d = 0; d < atoms.size(); ++d )
                            {
                                // итератор по контрактированным орбиталям четвертого атома
                                for ( size_t l = 0; l < atoms[d].get_element()->getCGOCount(); ++l )
                                {
                                    electronRepulsionTensor(it1, it2, it3, it4) = \
                     electronRepulsionCGO( atoms[a].get_element()->getCGO(i), atoms[a].getCoords(),
                                           atoms[b].get_element()->getCGO(j), atoms[b].getCoords(),
                                           atoms[c].get_element()->getCGO(k), atoms[c].getCoords(),
                                           atoms[d].get_element()->getCGO(l), atoms[d].getCoords() );
                                    ++it4;
                                }
                            }
                            ++it3;
                        }
                    }
                    ++it2;
                }
            }
            ++it1;
        }
    }
}

void Molecule::showElectronAttractionMatrix( )
{
    std::ofstream out("two-electron.txt");
    out << std::fixed << std::setprecision(8);

    for ( int i = 0; i < electronRepulsionTensor.dimension(0); ++i )
    {
        for ( int j = 0; j < electronRepulsionTensor.dimension(1); ++j )
        {
            for ( int k = 0; k < electronRepulsionTensor.dimension(2); ++k )
            {
                for ( int l = 0; l < electronRepulsionTensor.dimension(3); ++l )
                {
                    if ( std::abs(electronRepulsionTensor(i, j, k, l)) > 1e-10 )
                    {
                        out << "E[" << i << "," << j << "," << k << "," << l << "] = " << electronRepulsionTensor(i, j, k, l) << std::endl;
                    }
                }
            }
        }
    }

    out.close();
}


// приближение голых ядер, нулевая начальная матрица P
void Molecule::makeInitialGuess( )
{
    int size_ = size();
    matrixP = Eigen::MatrixXd::Zero( size_, size_ );
}

// Szabo, Ostlund.
// формула (3.154)
// G matrix -- two-electron part of the Fock matrix
void Molecule::fillGMatrix( )
{
    int size_ = size();

    // итерируемся по mu
    for ( int i = 0; i < size_; ++i )
    {
        // итерируемся по nu
        for ( int j = 0; j < size_; ++j )
        {
            matrixG(i, j) = 0.0;

            // итерируемся по lambda
            for ( int k = 0; k < size_; ++k )
            {
                // итерируемся по sigma
                for ( int l = 0; l < size_; ++l )
                {
                    // coulomb = (i j| k l)
                    // exchange = (i l| k j)
                    matrixG(i, j) += matrixP(k, l) * (electronRepulsionTensor(i, j, l, k) - 0.5 * electronRepulsionTensor(i, l, k, j));
                }
            }
        }
    }
}

void Molecule::SCF_initialize()
{
    int size_ = size();

    matrixF = Eigen::MatrixXd::Zero( size_, size_ );
    matrixFprime = Eigen::MatrixXd::Zero( size_, size_ );
    matrixG = Eigen::MatrixXd::Zero( size_, size_ );
    matrixC = Eigen::MatrixXd::Zero( size_, size_ );
    matrixP_new = Eigen::MatrixXd::Zero( size_, size_ );

    // filling Hcore matrix
    matrixHcore = kineticEnergyMatrix + nuclearAttractionMatrix;
    outFile << "Hcore matrix: " << std::endl << matrixHcore << std::endl << std::endl;
}

void Molecule::fillFMatrix()
{
    int size_ = size();
    for ( int i = 0; i < size_; ++i )
        for ( int j = 0; j < size_; ++j )
            matrixF(i, j) = matrixHcore(i, j) + matrixG(i, j);
}

void Molecule::setCharge(const int icharge)
// icharge - заряд молекулы
{
    int charge_ = 0;
    for ( size_t i = 0; i < atoms.size(); ++i )
        charge_ += atoms[i].get_element()->getCharge();

    charge = charge_ + icharge;
    set_charge = true;
}

void Molecule::SCF()
{
    if ( !set_charge )
        throw std::runtime_error("Molecule charge is not set.");

    std::vector<double> energy;

    int size_ = size();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;

    outFile << "-----------------------------------------------------" << std::endl;
    outFile << "Starting SCF procedure" << std::endl;

    for ( int it = 0; ; ++it )
    {
        if ( energy.size() > 2 )
            if ( std::abs(energy.end()[-1] - energy.end()[-2]) < 1e-12 )
                break;

        fillGMatrix();
        outFile << "matrix G: " << std::endl << matrixG << std::endl << std::endl;

        fillFMatrix();
        outFile << "Fock matrix: " << std::endl << matrixF << std::endl << std::endl;

        es.compute( matrixF, overlapMatrix );
        HF_OrbitalEnergies = es.eigenvalues();

        outFile << "Orbital energies: " << std::endl << HF_OrbitalEnergies << std::endl << std::endl;
        energy.push_back( HF_OrbitalEnergies(0) );

        matrixC = es.eigenvectors();

        outFile << "Coefficient matrix C: " << std::endl << matrixC << std::endl << std::endl;

        // mu
        for ( int i = 0; i < size_; ++i )
        {
            // nu
            for ( int j = 0; j < size_; ++j )
            {
                double res = 0.0;
                for ( int a = 0; a < charge / 2; ++a )
                    res += matrixC(i, a) * matrixC(j, a);

                matrixP_new(i, j) = 2 * res;
            }
        }

        outFile << "Density matrix P: " << std::endl << matrixP_new << std::endl << std::endl;
        matrixP = matrixP_new;

        // E_0; Szabo (3.184)
        // \mu
        double E_orbital = 0.0;
        double E_electronic = 0.0;
        for ( int i = 0; i < size_; ++i )
        {
            // \nu
            for ( int j = 0; j < size_; ++j )
            {
                E_orbital += 0.5 * matrixP(j, i) * matrixF(i, j);
                E_electronic += 0.5 * matrixP(j, i) * (kineticEnergyMatrix(i, j) + nuclearAttractionMatrix(i, j));
            }
        }
        outFile << "Orbital energy: " << E_orbital << std::endl;
        outFile << "Electronic energy: " << E_electronic << std::endl;

        double E0 = E_orbital + E_electronic;
        outFile << "E0: " << E0 << std::endl;

        // Enuc
        double Enuc = 0.0;
        for ( size_t i = 0; i < atoms.size(); ++i )
        {
            for ( size_t j = i + 1; j < atoms.size(); ++j )
            {
                Enuc += (atoms[i].get_element()->getCharge() * atoms[j].get_element()->getCharge()) / \
                    std::pow ( std::pow(atoms[i].get_x() - atoms[j].get_x(), 2) + \
                               std::pow(atoms[i].get_y() - atoms[j].get_y(), 2) + \
                               std::pow(atoms[i].get_z() - atoms[j].get_z(), 2), 0.5 );
            }
        }

        double Etot = E0 + Enuc;
        outFile << ">>> Etot: " << Etot << std::endl;

        outFile << "-------------------------------" << std::endl;
        outFile << "Iteration " << it << std::endl;
        outFile << "Highest orbital energy = " << energy.end()[-1] << std::endl;
        outFile << "Energy = " << Etot << std::endl;
        outFile << "-------------------------------" << std::endl << std::endl;
    }

	outFile << std::endl << "-------------------------------------------" << std::endl;
	outFile << "SCF procedure finished." << std::endl;
}

void Molecule::fillTwoElectronMOIntegrals()
{
    int size_ = size();
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
                    //std::cout << "twoElectronMOIntegrals(" << p << ", " << q << ", " << r << ", " << s << "): " <<
                    //             twoElectronMOIntegrals(p, q, r, s) << std::endl;
                }
            }
        }
    }
}

double Molecule::computeMP2_correction()
{
    double num = 0.0;
    double MP2_Energy = 0.0;

    int occ = charge / 2; // number of doubly occupied orbitals
    int size_ = size();

    std::cout << "(computeMP2_Energy) charge: " << charge << "; size_: " << size_ << std::endl;

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

void Molecule::fillAS_MO_TwoElectronIntegrals()
// Antisymmetrized integrals over spin-orbitals
{
    // количество спинорбиталей вдвое больше количества дважды занятых орбиталей
    int size_ = 2 * size();
    ASTwoElectronMOIntegrals.resize( size_, size_, size_, size_ );

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

void Molecule::fillSOHcoreMatrix()
{
    int size_ = 2 * size();
    SOHcoreMatrix.resize( size_, size_ );

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            SOHcoreMatrix(p, q) = matrixHcore(p/2, q/2) * (p%2 == q%2);
        }
    }
}

void Molecule::fillSOFockMatrix()
// fill spin-orbital Fock matrix
{
    int size_ = 2 * size();
    int occ = charge;
    SOFockMatrix.resize( size_, size_ );

    double res = 0;

    for ( int p = 0; p < size_; ++p )
    {
        for ( int q = 0; q < size_; ++q )
        {
            res = SOHcoreMatrix(p, q);
            for ( int m = 0; m < occ; ++m )
                res += ASTwoElectronMOIntegrals(p, m, q, m);

            SOFockMatrix(p, q) = res;
        }
    }
}

void Molecule::fill_initial_tia()
{
    int size_ = 2 * size();
    int occ = charge;
    tia.resize( size_, size_ );

    // i - occupied
    // a - unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int a = occ; a < size_; ++a )
        {
            tia(i, a) = 0.0;
        }
    }
}

void Molecule::fill_initial_tijab()
{
    int size_ = 2 * size();
    tijab.resize( size_, size_, size_, size_ );

    int occ = charge;

    // i, j -- occupied
    // a, b -- unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    // HF_OrbitalEnergies(i/2) -- потому что дважды занятые орбитали перешли в две спин-орбитали
                    tijab(i, j, a, b) = ASTwoElectronMOIntegrals(i, j, a, b) / \
                            (HF_OrbitalEnergies(i/2) + HF_OrbitalEnergies(j/2) - HF_OrbitalEnergies(a/2) - HF_OrbitalEnergies(b/2));
                }
            }
        }
    }
}

double Molecule::testCCSD_MP2_Energy()
{
    double res = 0.0;

    int size_ = 2 * size();
    int occ = charge;

    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    res += ASTwoElectronMOIntegrals(i, j, a, b) * tijab(i, j, a, b);
                    //std::cout << "twoelectronMOIN: " << ASTwoElectronMOIntegrals(i, j, a, b) << std::endl;
                    //std::cout << "t2: " << t2(i, j, a, b) << std::endl;
                }
            }
        }
    }

    return res / 4.0;
}

void Molecule::fill_tau_ijab()
{
    // i, j -- occupied
    // a, b -- unoccupied

    int size_ = size();
    int occ = charge;

    tau_ijab.resize( size_, size_, size_, size_ );

    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for( int b = occ; a < size_; ++b )
                {
                    tau_ijab(i, j, a, b) = tijab(i, j, a, b) + tia(i, a) * tia(j, b) -  tia(i, b) * tia(j, a);
                }
            }
        }
    }
}

void Molecule::fill_tau_tilda_ijab()
{
    int size_ = size();
    int occ = charge;

    tau_tilda_ijab.resize( size_, size_, size_, size_ );
    
    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    tau_tilda_ijab(i, j, a, b) = tijab(i, j, a, b) + 0.5 * (tia(i, a) * tia(j, b) -  tia(i, b) * tia(j, a));
                }
            }
        }
    }
}

void Molecule::fill_F()
{
    int size_ = 2 * size();
    int occ = charge;

    intermF.resize( size_, size_ );

   	// a, e -- unoccupied
    for ( int a = occ; a < size_; ++a )
    {
        for ( int e = occ; e < size_; ++e )
        {
            double sum1 = 0.0;
            // m -- occupied
            for ( int m = 0; m < occ; ++m )
            {
                sum1 += SOFockMatrix(m, e) * tia(m, a);
            }
            sum1 *= 0.5;

            double sum2 = 0.0;
            // m -- occupied
            // f -- unoccupied
            for ( int m = 0; m < occ; ++m )
            {
                for ( int f = occ; f < size_; ++f )
                {
                    sum2 += tia(f, m) * ASTwoElectronMOIntegrals(m, a, f, e);
                }
            }

            double sum3 = 0.0;
            // m,n -- occupied
            // f -- unoccupied
            for ( int m = 0; m < occ; ++m )
            {
                for ( int n = 0; n < occ; ++n )
                {
                    for ( int f = occ; f < size_; ++f )
                    {
                        sum3 += tau_tilda_ijab(m, n, a, f) * ASTwoElectronMOIntegrals(m, n, e, f);
                    }
                }
            }
            sum3 *= 0.5;

            intermF(a, e) = (1-delta(a, e)) * SOFockMatrix(a, e) - sum1 + sum2 - sum3;
        }
    }

    // i, m -- occupied
    for ( int m = 0; m < occ; ++m )
    {
        for ( int i = 0; i < occ; ++i )
        {
            double sum1 = 0.0;
            // e -- unoccupied
            for ( int e = occ; e < size_; ++e )
            {
                sum1 += tia(i, e) * SOFockMatrix(m, e);
            }
            sum1 *= 0.5;

            double sum2 = 0.0;
            // e -- unoccupied
            // n -- occupied
            for ( int e = occ; e < size_; ++e )
            {
                for ( int n = 0; n < occ; ++n )
                {
                    sum2 += tia(n, e) * ASTwoElectronMOIntegrals(m, n, i, e);
                }
            }

            double sum3 = 0.0;
            // n -- occupied
            // e, f -- unoccupied
            for ( int e = occ; e < size_; ++e )
            {
                for ( int f = occ; f < size_; ++f )
                {
                    for ( int n = 0; n < occ; ++n )
                    {
                        sum3 += tau_tilda_ijab(i, n, e, f) * ASTwoElectronMOIntegrals(m, n, e, f);
                    }
                }
            }
            sum3 *= 0.5;

            intermF(m, i) = (1-delta(m, i)) * SOFockMatrix(m, i) + sum1 + sum2 + sum3;
        }
    }

    // m -- occupied
    // e -- unoccupied
    for ( int m = 0; m < occ; ++m )
    {
        for ( int e = occ; e < size_; ++e )
        {
            double sum1 = 0.0;
            // n -- occupied
            // f -- unoccupied
            for ( int n = 0; n < occ; ++n )
            {
                for ( int f = occ; f < size_; ++f )
                {
                    sum1 += tia(n, f) * ASTwoElectronMOIntegrals(m, n, e, f);
                }
            }

            intermF(m, e) = SOFockMatrix(m, e) + sum1;
        }
    }
}

void Molecule::fill_W()
{
    int size_ = 2 * size();
    int occ = charge;

    intermW.resize( size_, size_, size_, size_ );

    // i, j, m, n -- occupied
    for ( int m = 0; m < occ; ++m )
    {
        for ( int n = 0; n < occ; ++n )
        {
            for ( int i = 0; i < occ; ++i )
            {
                for ( int j = 0; j < occ; ++j )
                {
                    double sum1 = 0.0;
                    // e -- unoccupied
                    for ( int e = occ; e < size_; ++e )
                    {
                        sum1 += tia(j, e) * ASTwoElectronMOIntegrals(m, n, i, e);
                    }

                    double sum2 = 0.0;
                    // e -- unoccupied
                    for ( int e = occ; e < size_; ++e )
                    {
                        sum2 += tia(i, e) * ASTwoElectronMOIntegrals(m, n, j, e);
                    }

                    double sum3 = 0.0;
                    // e, f -- unoccupied
                    for ( int e = occ; e < size_; ++e )
                    {
                        for ( int f = occ; f < size_; ++f )
                        {
                            sum3 += tau_ijab(i, j, e, f) * ASTwoElectronMOIntegrals(m, n, e, f);
                        }
                    }
                    sum3 *= 0.25;

                    intermW(m, n, i, j) = ASTwoElectronMOIntegrals(m, n, i, j) + sum1 - sum2 + sum3;
                }
            }
        }
    }

    // a, b, e, f -- unoccupied
    for ( int a = occ; a < size_; ++a )
    {
        for ( int b = occ; b < size_; ++b )
        {
            for ( int e = occ; e < size_; ++e )
            {
                for ( int f = occ; f < size_; ++f )
                {
                    double sum1 = 0.0;
                    // m -- occupied
                    for ( int m = 0; m < occ; ++m )
                    {
                        sum1 += tia(m, b) * ASTwoElectronMOIntegrals(a, m, e, f);
                    }

                    double sum2 = 0.0;
                    // m -- occupied
                    for ( int m = 0; m < occ; ++m )
                    {
                        sum2 += tia(m, a) * ASTwoElectronMOIntegrals(b, m, e, f);
                    }

                    double sum3 = 0.0;
                    // m, n -- occupied
                    for ( int m = 0; m < occ; ++m )
                    {
                        for ( int n = 0; n < occ; ++n )
                        {
                            sum3 += tau_ijab(m, n, a, b) * ASTwoElectronMOIntegrals(m, n, e, f);
                        }
                    }
                    sum3 *= 0.25;

                    intermW(a, b, e, f) = ASTwoElectronMOIntegrals(a, b, e, f) - sum1 + sum2 + sum3;
                }
            }
        }
    }

    // m, j -- occupied
    // b, e -- unccupied
    for ( int m = 0; m < occ; ++m )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int b = occ; b < size_; ++b )
            {
                for ( int e = occ; e < size_; ++e )
                {
                    double sum1 = 0.0;
                    // f -- unoccupied
                    for ( int f = occ; f < size_; ++f )
                    {
                        sum1 += tia(j, f) * ASTwoElectronMOIntegrals(m, b, e, f);
                    }

                    double sum2 = 0.0;
                    // n -- occupied
                    for ( int n = 0; n < occ; ++n )
                    {
                        sum2 += tia(n, b) * ASTwoElectronMOIntegrals(m, n, e, j);
                    }

                    double sum3 = 0.0;
                    // n -- occupied
                    // f -- unoccupied
                    for ( int n = 0; n < occ; ++n )
                    {
                        for ( int f = occ; f < size_; ++f )
                        {
                            sum3 += (0.5 * tijab(j, n, f, b) + tia(j, f) * tia(n, b)) * \
                                    ASTwoElectronMOIntegrals(m, n, e, f);
                        }
                    }

                    intermW(m, b, e, j) = ASTwoElectronMOIntegrals(m, b, e, j) + sum1 - sum2 - sum3;
                }
            }
        }
    }
}

void Molecule::fill_D1()
{
    int size_ = 2 * size();
    int occ = charge;

    D1.resize( size_, size_ );
    // a -- unoccupied, i -- occupied
    for ( int a = occ; a < size_; ++a )
    {
        for ( int i = 0; i < occ; ++i )
        {
            D1(i, a) = SOFockMatrix(i, i) - SOFockMatrix(a, a);
        }
    }
}

void Molecule::fill_D2()
{
    int size_ = 2 * size();
    int occ = charge;

    D2.resize( size_, size_, size_, size_ );
    // i, j -- occupied
    // a, b -- unoccupied

    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    D2(i, j, a, b) = SOFockMatrix(i, i) + SOFockMatrix(j, j) - SOFockMatrix(a, a) - SOFockMatrix(b, b);
                }
            }
        }
    }
}

void Molecule::update_t1()
{
    int size_ = 2 * size();
    int occ = charge;

    t1_updated.resize( size_, size_ );

    // i -- occupied
    // a -- unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int a = occ; a < size_; ++a )
        {
            double sum1 = 0.0;
            // e -- unoccupied
            for ( int e = occ; e < size_; ++e )
            {
                sum1 += tia(i, e) * intermF(a, e);
            }

            double sum2 = 0.0;
            // m -- occupied
            for ( int m = 0; m < occ; ++m )
            {
                sum2 += tia(m, a) * intermF(m, i);
            }

            double sum3 = 0.0;
            // m -- occupied, e -- unoccupied
            for ( int m = 0; m < occ; ++m )
            {
                for ( int e = occ; e < size_; ++e )
                {
                    sum3 += tijab(i, m, a, e) * intermF(m, e);
                }
            }

            double sum4 = 0.0;
            // n -- occupied, f -- unoccupied
            for ( int n = 0; n < occ; ++n )
            {
                for ( int f = occ; f < size_; ++f )
                {
                    sum4 += tia(n, f) * ASTwoElectronMOIntegrals(n, a, i, f);
                }
            }

            double sum5 = 0.0;
            // m -- occupied; e, f -- unoccupied
            for ( int m = 0; m < occ; ++m )
            {
                for ( int e = occ; e < size_; ++e )
                {
                    for ( int f = occ; f < size_; ++f )
                    {
                        sum5 += tijab(i, m, e, f) * ASTwoElectronMOIntegrals(m, a, e, f);
                    }
                }
            }
            sum5 *= 0.5;

            double sum6 = 0.0;
            // m, n -- occupied; e -- unoccupied
            for ( int m = 0; m < occ; ++m )
            {
                for ( int n = 0; n < occ; ++n )
                {
                    for ( int e = occ; e < size_; ++e )
                    {
                        sum6 += tijab(m, n, a, e) * ASTwoElectronMOIntegrals(n, m, e, i);
                    }
                }
            }
            sum6 *= 0.5;

            t1_updated(i, a) = SOFockMatrix(i, a) + sum1 - sum2 + sum3 - sum4 - sum5 - sum6;
            t1_updated(i, a) /= D1(i, a);
        }
    }
}

void Molecule::update_t2()
{
    int size_ = 2 * size();
    int occ = charge;

    t2_updated.resize( size_, size_, size_, size_ );

    // i, j -- occupied; a, b -- unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    double sum1ab = 0.0;
                    // меняем местами a и b (действие оператора P_(a, b)
                    double sum1ba = 0.0;

                    // e -- unoccupied
                    for ( int e = occ; e < size_; ++e )
                    {
                        double res = 0.0;
                        // m -- occupied
                        for ( int m = 0; m < occ; ++m )
                        {
                            res += tia(m, b) * intermF(m, e);
                        }
                        sum1ab += tijab(i, j, a, e) * (intermF(b, e) - 0.5 * res);

                        res = 0.0;
                        // m -- occupied
                        for ( int m = 0; m < occ; ++m )
                        {
                            res += tia(m, a) * intermF(m, e);
                        }
                        sum1ba += tijab(i, j, b, e) * (intermF(a, e) - 0.5 * res);
                    }

                    double sum2ij = 0.0;
                    // меняем местами i и j (действие оператора P_(i, j)
                    double sum2ji = 0.0;

                    // m -- occupied
                    for ( int m = 0; m < occ; ++m )
                    {
                        double res = 0.0;
                        // e -- unoccupied
                        for ( int e = occ; e < size_; ++e )
                        {
                            res += tia(j, e) * intermF(m, e);
                        }

                        sum2ij += tijab(i, m, a, b) * (intermF(m, j) + 0.5 * res);

                        res = 0.0;
                        for ( int e = occ; e < size_; ++e )
                        {
                            res += tia(i, e) * intermF(m, e);
                        }

                        sum2ji += tijab(j, m, a, b) * (intermF(m, i) + 0.5 * res);
                    }

                    double sum3 = 0.0;
                    // m, n -- occupied
                    for ( int m = 0; m < occ; ++m )
                    {
                        for ( int n = 0; n < occ; ++n )
                        {
                            sum3 += tau_ijab(m, n, a, b) * intermW(m, n, i, j);
                        }
                    }
                    sum3 *= 0.5;

                    double sum4 = 0.0;
                    // e, f -- unoccupied
                    for ( int e = occ; e < size_; ++e )
                    {
                        for ( int f = occ; f < size_; ++f )
                        {
                            sum4 += tau_ijab(i, j, e, f) * intermW(a, b, e, f);
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

                    for ( int m = 0; m < occ; ++m )
                    {
                        for ( int e = occ; e < size_; ++e )
                        {
                            sum5ijab += (tijab(i, m, a, e) * intermW(m, b, e, j) - \
                                         tia(i, e) * tia(m, a) * ASTwoElectronMOIntegrals(m, b, e, j));
                            sum5ijba += (tijab(i, m, b, e) * intermW(m, a, e, j) - \
                                         tia(i, e) * tia(m, b) * ASTwoElectronMOIntegrals(m, a, e, j));
                            sum5jiab += (tijab(j, m, a, e) * intermW(m, b, e, i) - \
                                         tia(j, e) * tia(m, a) * ASTwoElectronMOIntegrals(m, b, e, i));
                            sum5jiba += (tijab(j, m, b, e) * intermW(m, a, e, i) - \
                                         tia(j, e) * tia(m, b) * ASTwoElectronMOIntegrals(m, a, e, i));
                        }
                    }

                    // 2 суммы, оператор P_(i, j)
                    double sum6ij = 0.0;
                    // переставляем местами (i, j)
                    double sum6ji = 0.0;

                    // e -- unoccupied
                    for ( int e = occ; e < size_; ++e )
                    {
                        sum6ij += (tia(i, e) * ASTwoElectronMOIntegrals(a, b, e, j));
                        sum6ji += (tia(j, e) * ASTwoElectronMOIntegrals(a, b, e, i));
                    }

                    // 2 суммы, оператор P_(a, b)
                    double sum7ab = 0.0;
                    // переставляем местами (a, b)
                    double sum7ba = 0.0;

                    // m -- occupied
                    for ( int m = 0; m < occ; ++m )
                    {
                        sum7ab += (tia(m, a) * ASTwoElectronMOIntegrals(m, b, i, j));
                        sum7ba += (tia(m, b) * ASTwoElectronMOIntegrals(m, a, i, j));
                    }

                    t2_updated(i, j, a, b) = ASTwoElectronMOIntegrals(i, j, a, b) + sum1ab - sum1ba - sum2ij + sum2ji + \
                        sum3 + sum4 + sum5ijab - sum5ijba - sum5jiab + sum5jiba + sum6ij - sum6ji - sum7ab + sum7ba;
                    t2_updated(i, j, a, b) /= D2(i, j, a, b);
                }
            }
        }
    }
}

double Molecule::computeCCSD_correction()
{
    int size_ = 2 * size();
    int occ = charge;

    double sum1 = 0.0;
    // i -- occupied, a -- unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int a = occ; a < size_; ++a )
        {
            sum1 += SOFockMatrix(i, a) * t1_updated(i, a);
        }
    }
    std::cout << "(CCSD correction) sum1: " << sum1 << std::endl;

    double sum2 = 0.0;
    // i, j -- occupied; a, b -- unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    sum2 += ASTwoElectronMOIntegrals(i, j, a, b) * t2_updated(i, j, a, b);
                }
            }
        }
    }
    sum2 *= 0.25;
    std::cout << "(CCSD correction) sum2: " << sum2 << std::endl;

    double sum3 = 0.0;
    // i, j -- occupied; a, b -- unoccupied
    for ( int i = 0; i < occ; ++i )
    {
        for ( int j = 0; j < occ; ++j )
        {
            for ( int a = occ; a < size_; ++a )
            {
                for ( int b = occ; b < size_; ++b )
                {
                    sum3 += ASTwoElectronMOIntegrals(i, j, a, b) * t1_updated(i, a) * t1_updated(j, b);
                }
            }
        }
    }
    sum3 *= 0.5;
    std::cout << "(CCSD correction) sum3: " << sum3 << std::endl;

    return sum1 + sum2 + sum3;
}
