#include "Element.hpp"

void Element::add_basis_functions( BasisFunctionBuilder * bfb )
{
	for ( auto * bf : bfb->getBasisFunctions() )
		basisFunctions.push_back( bf );
}


void Element::show()
{
	std::cout << "Element name: " << name << std::endl;
	for ( auto bf : basisFunctions )
		bf->show();
}

