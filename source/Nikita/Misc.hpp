#pragma once

#include <cmath>     // библиотека математических функций
#include <boost/math/special_functions/gamma.hpp>
#include "CoordsClass.hpp"

int d_factorial( int n ){
	return ( n <= 1 ) ? 1 : n * d_factorial( n - 2 );
}

double factor( int i, double alpha ){
	return  d_factorial( 2 * i - 1 ) * std::sqrt( M_PI ) / std::pow( 2, i ) / std::pow( alpha, i + 0.5 ) ;
}

double boys_function( const int n, const double x ){
	if ( std::abs(x) < 10e-10 ) { return 1.0 / ( 2.0 * n + 1.0 ); }
	else { 	return 0.5 * std::pow(x,-0.5-n) * ( std::tgamma( 0.5 + n ) - boost::math::tgamma( 0.5 + n, x ) ); }
}

_Coords gauss_prod_center( double a, _Coords A,\
			   double b, _Coords B ){
	double p = 1 / (a + b);
	double x = ( a * A.get_x() + b * B.get_x() ) * p;
	double y = ( a * A.get_y() + b * B.get_y() ) * p;
	double z = ( a * A.get_z() + b * B.get_z() ) * p;
	return _Coords(x,y,z);
}

bool unequal( _Coords c1, _Coords c2 ){
	if ( ( c1.get_x() != c2.get_x() ) || \
	     ( c1.get_y() != c2.get_y() ) || \
	     ( c1.get_z() != c2.get_z() ) ){
		return true;
	}
}
