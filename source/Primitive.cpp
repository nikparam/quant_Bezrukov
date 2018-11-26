#include "Primitive.hpp"

void Primitive::normalize( QuantumNumbers const & qNumbers )
{
	int L = qNumbers.getAngularMomentum();
    double norm = std::pow( 2.0, 2 * L + 1.5 ) * std::pow( exponent, L + 1.5 ) / MathUtils::doubleFactorial(2 * qNumbers.i - 1) / MathUtils::doubleFactorial(2 * qNumbers.j - 1) / MathUtils::doubleFactorial(2 * qNumbers.k - 1) / std::pow(M_PI, 1.5);
	coefficient *= std::pow(norm, 0.5);	
}

double Primitive::get_norm( QuantumNumbers const& qNumbers ) const
{
    int L = qNumbers.getAngularMomentum();
    return std::pow( 2.0, 2 * L + 1.5 ) * std::pow( exponent, L + 1.5 ) / MathUtils::doubleFactorial(2 * qNumbers.i - 1) / MathUtils::doubleFactorial(2 * qNumbers.j - 1) / MathUtils::doubleFactorial(2 * qNumbers.k - 1) / std::pow(M_PI, 1.5);
}

void Primitive::show()
{
	std::cout << "(Primitive) exponent = " << exponent << "; coefficient = " << coefficient << std::endl; 
}
