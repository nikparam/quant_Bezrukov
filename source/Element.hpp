#pragma once

#include <vector>
#include <string>
#include "BasisFunctionBuilder.hpp"

class Element
{
public:
	Element( std::string name ) : name(name)
	{
	}

	~Element()
	{
		for ( BasisFunction * bf : basisFunctions )
			delete bf;
	}

	void add_basis_functions( BasisFunctionBuilder * bfb );

	void show();

	BasisFunction * getBasisFunction( const int n ) { return basisFunctions.at( n ); }

	size_t getBasisFunctionsCount() { return basisFunctions.size(); }

private:
	std::string name;
	std::vector<BasisFunction*> basisFunctions;
};

