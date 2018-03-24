#pragma once

#include <cmath>     // библиотека математических функций

int d_factorial( int n ){
	return ( n <= 1 ) ? 1 : n * d_factorial( n - 2 );
}

double factor( int i, double alpha ){
	return  d_factorial( 2 * i - 1 ) * std::sqrt( M_PI ) / std::pow( 2, i ) / std::pow( alpha, i + 0.5 ) ;
}

