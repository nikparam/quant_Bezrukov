#include "ContractedGaussianOrbital.hpp"

void ContractedGaussianOrbital::add_primitive( double alpha, double exponent )
{
	primitives.emplace_back( alpha, exponent );
	primitives.end()[-1].normalize( qNumbers );
}


void ContractedGaussianOrbital::show()
{
    std::cout << "(CGO)" << primitives.size() << std::endl;
    std::cout << "(CGO) QuantumNumbers: " << qNumbers.i << " " << qNumbers.j << " " << qNumbers.k << std::endl;

	for ( auto p : primitives )
		p.show();
}
