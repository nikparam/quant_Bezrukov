#pragma once

#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

class BasisFunction
{
public:
	BasisFunction( char angular_part )
		: angular_part(angular_part)
	{
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

	void showQuantumNumbers()
	{
		cout << "angular part: " << angular_part << endl;
		for ( auto v : quantumNumbers  )
			cout << v.i << " " << v.j << " " << v.k << endl;
		cout << endl << endl;
	}

	void add_primitive( int n, double exponent, double coefficient )
	{
		primitives.emplace_back( n, exponent, coefficient );
		primitives.end()[-1].renormalize();
	}

	Primitive * last_primitive ( )
	{
		return &primitives.end()[-1];
	}	

	void show()
	{
		cout << "Angular part: " << angular_part << endl;
		cout << "Number of added primitives: " << primitives.size() << endl;
	}

private:
	char angular_part;
	vector <Primitive> primitives;
	vector <QuantumNumbers> quantumNumbers;
};

