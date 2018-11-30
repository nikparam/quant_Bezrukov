#pragma once

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

extern "C" void boys_func_(int&, double&, double&);

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

    inline double BoysFunction( int n, double x )
    // ERROR
    {
        //if ( std::abs(x) < 1e-10 )
        //    return 1.0 / (2.0 * n + 1.0);
        //return 0.5 * std::pow(x, -0.5 - n) * (std::tgamma(0.5 + n) - boost::math::tgamma(0.5 + n, x));
        double res;
        boys_func_(n, x, res);

        return res;
    }
}
