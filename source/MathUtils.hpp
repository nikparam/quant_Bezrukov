#pragma once

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

namespace MathUtils
{
	inline unsigned int doubleFactorial( unsigned int n )
	{
		if ( n <=  1 ) return 1;

		unsigned int res = 1;
		for ( int k = n; k > 0; k -= 2 )
			res *= k;

		return res;
    }

    inline double BoysFunction( const int n, const double x )
    {
        if ( std::abs(x) < 1e-10 )
            return 1.0 / (2.0 * n + 1.0);
        return 0.5 * std::pow(x, -0.5 - n) * (std::tgamma(0.5 + n) - boost::math::tgamma(0.5 + n, x));
    }
}
