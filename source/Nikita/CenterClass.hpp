#pragma once

#include <iostream>  // библиотека стандартного ввода/вывода
#include <vector>    // библиотека для класса векторов

class _Center{

public:
	_Center( double x, double y, double z ){
		coords.push_back( x );
		coords.push_back( y );
		coords.push_back( z );
	}
private:
	std::vector<double> coords;

};

