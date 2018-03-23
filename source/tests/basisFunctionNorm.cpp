#include "basisFunctionNorm.hpp"

double checkNorm( BasisFunction * bf )
{
	double res = 0;

	unsigned int i = bf->getQuantumNumbers().i;
	unsigned int j = bf->getQuantumNumbers().j;
	unsigned int k = bf->getQuantumNumbers().k;

	unsigned int L = i + j + k;

	for ( size_t alpha = 0; alpha < bf->getPrimitivesCount(); alpha++ )
	{
		for ( size_t beta = alpha; beta < bf->getPrimitivesCount(); beta++ )
		{
			{
				double temp =  bf->getPrimitive(alpha).get_coefficient() * bf->getPrimitive(beta).get_coefficient() * pow(M_PI, 1.5) / pow(2.0, L) / pow(bf->getPrimitive(alpha).get_exponent() + bf->getPrimitive(beta).get_exponent(), L + 1.5);
			   	double fact = MathUtils::doubleFactorial(2 * i - 1) * MathUtils::doubleFactorial(2 * j - 1) * MathUtils::doubleFactorial(2 * k - 1);	  

				// диагональные члены учитываем один раз
				// а внедиагональные нужно прибавить дважды
				if ( alpha == beta )
					res += (temp * fact);
				else
					res += (2 * temp * fact);
			}
		}
	}

	return res;
}

