#pragma once

#include <iostream>  // библиотека стандартного ввода/вывода
#include <vector>    // библиотека для класса векторов

class _Coords{

public:
	_Coords( double x, double y, double z ){
		 coords.push_back( x );
		 coords.push_back( y );
		 coords.push_back( z );
	}

	double get_x() { return coords[0]; }
	double get_y() { return coords[1]; }
	double get_z() { return coords[2]; }

private:
	std::vector<double> coords;

};

