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
    DEBUG_METHOD("calculateEijt");
    std::cout << "i: " << i << "; j: " << j << "; t: " << t << "; Qx: " << Qx << "; a: " << a << "; b: " << b << std::endl;

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
    DEBUG_METHOD("electronRepulsionPrimitive");
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
    DEBUG_METHOD("electronRepulsionCGO");

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
    outFile << std::fixed << std::setprecision(12);
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
    //DEBUG_METHOD("fillElectronRepulsionMatrix");

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

void Molecule::showElectronRepulsionTensor( )
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
                  
                    // Какого черта здесь стояли коэффициенты 1 и 0.5, а не 2 и 1??
                    matrixG(i, j) += matrixP(k, l) * (2.0 * electronRepulsionTensor(i, j, l, k) - electronRepulsionTensor(i, l, k, j));
                    //matrixG(i, j) += matrixP(k, l) * (electronRepulsionTensor(i, j, l, k) - 0.5 * electronRepulsionTensor(i, l, k, j));
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
    //outFile << "Hcore matrix: " << std::endl << matrixHcore << std::endl << std::endl;
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

double Molecule::find_max( const Eigen::MatrixXd & B )
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

double Molecule::compute_Enuc()
{
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
    
    return Enuc;
}

double Molecule::compute_E0()
{
    int size_ = size();

    // E_0; Szabo (3.184)
    // \mu
    double E_orbital = 0.0;
    double E_electronic = 0.0;
    for ( int i = 0; i < size_; ++i )
    {
        // \nu
        for ( int j = 0; j < size_; ++j )
        {
            E_orbital += matrixP(j, i) * matrixF(i, j);
            E_electronic += matrixP(j, i) * (kineticEnergyMatrix(i, j) + nuclearAttractionMatrix(i, j));
        }
    }

    return E_orbital + E_electronic;
}

double Molecule::SCF_DIIS()
{
    if ( !set_charge )
        throw std::runtime_error("Molecule charge is not set.");

    int size_ = size();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;

    double Enuc = compute_Enuc();
    double E0 = 0.0;
    double Etot = 0.0;

    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << "Starting SCF procedure" << std::endl;
    
    fillGMatrix();
    fillFMatrix();
    es.compute( matrixF, overlapMatrix );
    matrixC = es.eigenvectors();
    for ( int i = 0; i < size_; ++i )
    {
        // nu
        for ( int j = 0; j < size_; ++j )
        {
            double res = 0.0;
            for ( int a = 0; a < charge / 2; ++a )
                res += matrixC(i, a) * matrixC(j, a);

            matrixP(i, j) = res; 
        }
    }

    Eigen::MatrixXd A = overlapMatrix.pow(-0.5);
    //std::cout << "A:\n" << A << std::endl;

    Eigen::MatrixXd diis_error;
    std::vector<Eigen::MatrixXd> FockMatrices;
    std::vector<Eigen::MatrixXd> DIIS_errors;

    double Etot_old = 0.0;
    for ( int it = 0; it < max_diis_iter; ++it, Etot_old = Etot )
    {
        fillGMatrix();
        fillFMatrix();

        //std::cout << "C:\n" << matrixC << std::endl;
        //std::cout << "D:\n" << matrixP << std::endl;
        //std::cout << "D (new):\n" << matrixP_new << std::endl;
        //std::cout << "S:\n" << overlapMatrix << std::endl;
        //std::cout << "F:\n" << matrixF << std::endl;

        diis_error = matrixF * matrixP * overlapMatrix - overlapMatrix * matrixP * matrixF;
        diis_error = A * diis_error * A; // wtf?
       
        // добавляем матрицу Фока и матрицу 'ошибок' в вектора
        FockMatrices.push_back( matrixF );
        DIIS_errors.push_back( diis_error );

        //std::cout << "Error matrix: \n" << diis_error << std::endl;

        E0 = compute_E0();
        Etot = E0 + Enuc;
        //std::cout << "E0: " << E0 << "; Enuc: " << Enuc << std::endl;
        std::cout << "SCF iteration: " << it << "; Energy: " << Etot << "; dE: " << Etot - Etot_old << std::endl; 
        
        if ( std::abs(Etot - Etot_old) < E_conv )
            break;

        if ( it >= 1 )
        {
            size_t diis_size = FockMatrices.size();
            //std::cout << "diis_size: " << diis_size << std::endl;
            // удаляем лишние, если diis_count > diis_max_size
            if ( diis_size > max_diis_size )
            {
                FockMatrices.erase( FockMatrices.begin() );
                DIIS_errors.erase( DIIS_errors.begin() );
                diis_size = max_diis_size;
            }
            //std::cout << "(after erasing) diis_size: " << diis_size << std::endl;

            // инициализируем матрицу DIIS минус единицами
            Eigen::MatrixXd B = Eigen::MatrixXd::Constant( diis_size + 1, diis_size + 1, -1.0 );
            B( diis_size, diis_size ) = 0.0;
          
            for ( size_t k = 0; k < diis_size; ++k )
            {
                for ( size_t l = 0; l < diis_size; ++l )
                {
                    double res = 0.0;
                    for ( int i = 0; i < DIIS_errors[k].rows(); ++i )
                        for ( int j = 0; j < DIIS_errors[k].cols(); ++j )
                            res += DIIS_errors[k](i, j) * DIIS_errors[l](i, j);

                    B(k, l) = res;
                }
            }

            // normalize 
            double max_ = find_max( B ); 
            for ( int i = 0; i < B.rows() - 1; ++i )
                for ( int j = 0; j < B.cols() - 1; ++j )
                    B(i, j) /= max_;

            //std::cout << "B:\n" << B << std::endl;
            
            Eigen::VectorXd resid = Eigen::VectorXd::Zero( diis_size + 1 );
            resid( diis_size ) = -1.0;

            Eigen::VectorXd ci = B.colPivHouseholderQr().solve( resid );
            //std::cout << "ci: \n" << ci << std::endl;

            //Eigen::MatrixXd F_new = Eigen::MatrixXd::Zero( matrixF.rows(), matrixF.cols() );
            matrixF = Eigen::MatrixXd::Zero( matrixF.rows(), matrixF.cols() );
            for ( size_t k = 0; k < diis_size; ++k )
            {
                //F_new += ci(k) * FockMatrices[k];
                matrixF += ci(k) * FockMatrices[k];
                //std::cout << "adding with coefficient " << ci(k) << std::endl << FockMatrices[k] << std::endl;
            }

            //std::cout << "F:\n" << matrixF << std::endl;
        }
        
        es.compute( matrixF, overlapMatrix );
        HF_OrbitalEnergies = es.eigenvalues();
        matrixC = es.eigenvectors();

        // Рассчитываем новую матрицу плотности
        // mu
        for ( int i = 0; i < size_; ++i )
        {
            // nu
            for ( int j = 0; j < size_; ++j )
            {
                double res = 0.0;
                for ( int a = 0; a < charge / 2; ++a )
                    res += matrixC(i, a) * matrixC(j, a);

                matrixP(i, j) = res; //2 * res;
            }
        }
    }

    return Etot;
}

double Molecule::SCF()
{
    if ( !set_charge )
        throw std::runtime_error("Molecule charge is not set.");

    int size_ = size();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;

    double E0 = 0.0;
    double Etot_old = 0.0;
    double Etot = 0.0;

    outFile << "-----------------------------------------------------" << std::endl;
    outFile << "Starting SCF procedure" << std::endl;

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

    for ( int it = 0; ; ++it, Etot_old = Etot )
    {
        fillGMatrix();
        fillFMatrix();
        //outFile << "Fock matrix: " << std::endl << matrixF << std::endl << std::endl;

        es.compute( matrixF, overlapMatrix );
        HF_OrbitalEnergies = es.eigenvalues();
        matrixC = es.eigenvectors();

        //outFile << "Coefficient matrix C: " << std::endl << matrixC << std::endl << std::endl;

        // mu
        for ( int i = 0; i < size_; ++i )
        {
            // nu
            for ( int j = 0; j < size_; ++j )
            {
                double res = 0.0;
                for ( int a = 0; a < charge / 2; ++a )
                    res += matrixC(i, a) * matrixC(j, a);

                matrixP_new(i, j) = res; //2 * res;
            }
        }

        //outFile << "Density matrix P: " << std::endl << matrixP_new << std::endl << std::endl;
        matrixP = matrixP_new;

        E0 = compute_E0();
        Etot = E0 + Enuc;

        if ( std::abs(Etot - Etot_old) < E_conv )
            break;

        std::cout << "SCF iteration: " << it << "; Energy: " << Etot << "; dE: " << Etot - Etot_old << std::endl;

        //outFile << "-------------------------------" << std::endl;
        //outFile << "Iteration " << it << std::endl;
        //outFile << "Highest orbital energy = " << energy.end()[-1] << std::endl;
        //outFile << "Energy = " << Etot << std::endl;
        //outFile << "-------------------------------" << std::endl << std::endl;
    }

    //outFile << std::endl << "-------------------------------------------" << std::endl;
    //outFile << "SCF procedure finished." << std::endl;

    return Etot;
}
