#include "BasisFunctionBuilder.hpp"

BasisFunctionBuilder::BasisFunctionBuilder( const char angularPart ) : angularPart(angularPart)	
{
	std::vector<QuantumNumbers> qNumbers = generateQuantumNumbers( angularPart );
	for ( auto & triple : qNumbers )
	{
		BasisFunction * bf = new BasisFunction( triple );
		basisFunctions.push_back( bf );
	}	
}
std::vector<QuantumNumbers> BasisFunctionBuilder::generateQuantumNumbers( const char angularPart )
{
	std::vector<QuantumNumbers> qNumbers;

	int l;

	switch( angularPart )
	{
		case 'S':
		{
			l = 0;
			break;
		}
		case 'P':
		{
			l = 1;
			break;	
		}
		case 'D':
		{
			l = 2;
			break;
		}
		case 'F':
		{
			l = 3;
			break;
		}
		default:
		{
			throw std::invalid_argument( "Angular part is not implemented" );
		}
	}
	
	for ( int i = 0; i <= l; i++ )
		for ( int j = 0; j <= l - i; j++ )
			for ( int k = 0; i + j + k <= l; k++ )
				if ( i + j + k == l )
					qNumbers.emplace_back( i, j, k );

	return qNumbers;
}	

void BasisFunctionBuilder::add_primitive( double alpha, double coeff )
{
	for ( auto * bf : basisFunctions )
		bf->add_primitive( alpha, coeff );
}
