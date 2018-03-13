#pragma once

#include <vector>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::vector;

class BasisFunction
{
public:
	BasisFunction( char angular_part )
		: angular_part(angular_part)
	{
		fillQuantumNumbers();
	}

	~BasisFunction()
	{
	}
	
	void fillQuantumNumbers()
	{
		int l;

		switch( angular_part )
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
						quantumNumbers.emplace_back( i, j, k );
	}

	void normalize( )
	{
		double sum = 0;

		for ( size_t i = 0; i < primitives.size(); i++ )
		{
			for ( size_t  j = i; j < primitives.size(); j++ )
			{
				if ( (primitives[i].get_qNumbers().i + primitives[j].get_qNumbers().i) % 2 == 0 &&
					 (primitives[i].get_qNumbers().j + primitives[j].get_qNumbers().j) % 2 == 0 &&
					 (primitives[i].get_qNumbers().k + primitives[j].get_qNumbers().k) % 2 == 0 )
				{
						//cout << "(first primitive " << i << ") " <<  primitives[i].get_qNumbers().i << " " << primitives[i].get_qNumbers().j << " " << primitives[i].get_qNumbers().k << endl;
						//cout << "(second primitive " << j << ") " << primitives[j].get_qNumbers().i << " " << primitives[j].get_qNumbers().j << " " << primitives[j].get_qNumbers().k << endl;
					double r = 0.5 * (primitives[i].get_qNumbers().getAngularMomentum() + primitives[j].get_qNumbers().getAngularMomentum());
					double temp =  primitives[i].get_coefficient() * primitives[j].get_coefficient() * pow(M_PI, 1.5) / pow(2.0, r) / pow(primitives[i].get_exponent() + primitives[j].get_exponent(), r + 1.5);
				   	double fact = MathUtils::doubleFactorial(primitives[i].get_qNumbers().i + primitives[j].get_qNumbers().i - 1) * 
						 		  MathUtils::doubleFactorial(primitives[i].get_qNumbers().j + primitives[j].get_qNumbers().j - 1) *
							  	  MathUtils::doubleFactorial(primitives[i].get_qNumbers().k + primitives[j].get_qNumbers().k - 1);	  
					//cout << "overlap: " << pow(M_PI, 1.5) / pow(2.0, r) / pow(primitives[i].get_exponent() + primitives[j].get_exponent(), r + 1.5) << endl;
					//cout << "temp: " << temp << endl;

					// диагональные члены учитываем один раз
					// а внедиагональные нужно прибавить дважды
					if ( i == j )
						sum = sum + temp * fact;
					else
						sum = sum + 2 * temp * fact;
				}
			}
		}

		cout << std::fixed << std::setprecision(5);
		cout << "(BasisFunction) sum = " << sum << endl;
	}

	void showQuantumNumbers()
	{
		cout << "angular part: " << angular_part << endl;
		for ( auto v : quantumNumbers  )
			cout << v.i << " " << v.j << " " << v.k << endl;
		cout << endl << endl;
	}

	void add_primitive( double exponent, double coefficient )
	{
		// добавляется примитив, тут же перенормируется, т.к. для этого есть все необходимое 
		for ( auto & qN : quantumNumbers )
			primitives.emplace_back( exponent, coefficient, qN );
	}

	void show()
	{
		cout << "Angular part: " << angular_part << endl;
		cout << "Number of added primitives: " << primitives.size() << endl;
		for ( auto p : primitives )
			p.show();

		normalize();
	}

private:
	char angular_part;
	vector <Primitive> primitives;
	vector <QuantumNumbers> quantumNumbers;
};

