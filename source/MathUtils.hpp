#pragma once

namespace MathUtils
{
	unsigned int doubleFactorial( unsigned int n )
	{
		if ( n <=  1 ) return 1;

		unsigned int res = 1;
		for ( int k = n; k > 0; k -= 2 )
			res *= k;

		return res;
	}
}
