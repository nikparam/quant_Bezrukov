#include "../QuantumNumbers.hpp"
#include "../Primitive.hpp"
#include "../BasisFunction.hpp"
#include "../Element.hpp"
#include "../Basis.hpp"

void checkNorm( BasisFunction * bf )
{
	double sum = 0;

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
					sum = sum + temp * fact;
				else
					sum = sum + 2 * temp * fact;
			}
		}
	}

	cout << "(checkNorm) scalar product = " << sum << endl;
}


int main()
{
	Basis basis;
	basis.read("../basis/cc-pvdz.gamess-us.dat");

	std::cout << std::fixed << std::setprecision(8);

	for ( size_t i = 0; i < basis.getElementsCount(); i++ )
	{
		for ( size_t j = 0; j < basis.getElement(i)->getBasisFunctionsCount(); j++ )
		{
			checkNorm( basis.getElement(i)->getBasisFunction(j) );
		}
	}

	return 0;
}

