#pragma once
#include <iostream>
#include <cmath>

#include "CoordsClass.hpp"

class _Atom{

public:
	_Atom( double x, double y, double z, _Element * e): coords(x,y,z), e(e) { }
	~_Atom(){}

	void show_atom(){
		std::cout << e -> get_name() << ": " \
			  << ( coords.get_x() ) << " " \
			  << ( coords.get_y() ) << " " \
			  << ( coords.get_z() ) << std::endl;
		e -> show_basis();
	}

	_Coords get_c(){ return coords; }
	_Element * get_e(){ return e; }

private:
	_Coords coords;
	_Element * e;

};
