#include "molecule.h"

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
    //std::cout << "a1: " << a1->get_name() << "; " << a1->get_x() << " " << a1->get_y() << " " << a1->get_z() << std::endl;
    //std::cout << "a2: " << a2->get_name() << "; " << a2->get_x() << " " << a2->get_y() << " " << a2->get_z() << std::endl;

    size_t size = 0;
    for ( size_t i = 0; i < atoms.size(); i++ )
        size += atoms[i].get_element()->getCGOCount();
    // std::cout << "Size of overlap matrix: " << size << " x " << size << std::endl;

    overlapMatrix.resize(size, size);

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
                    if ( it2 == size )
                    {
                        ++it1;
                        it2 = 0;
                    }

                }
            }
        }
    }
}

void Molecule::saveOverlapMatrix( std::string filename )
{
    if ( filename == "" )
        filename = "__overlapMatrix.dat";

    std::ofstream outFile( filename );
    outFile << std::fixed << std::setprecision(8);
    outFile << overlapMatrix << std::endl;

    outFile.close();
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
