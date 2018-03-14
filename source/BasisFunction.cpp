#include "BasisFunction.hpp"

void BasisFunction::add_primitive( double alpha, double exponent )
{
	primitives.emplace_back( alpha, exponent );
	primitives.end()[-1].normalize( qNumbers );
}


void BasisFunction::show()
{
	std::cout << "(BasisFunction) Number of primitives: " << primitives.size() << std::endl;
	std::cout << "(BasisFunction) QuantumNumbers: " << qNumbers.i << " " << qNumbers.j << " " << qNumbers.k << std::endl;

	for ( auto p : primitives )
		p.show();
}
