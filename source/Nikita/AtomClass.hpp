#pragma once
#include <iostream>
#include <cmath>

#include "CoordsClass.hpp"

class _Atom{

public:
	_Atom( double x, double y, double z, _Element * e): coords(x,y,z), e(e) { add_charge(); }
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

	void add_charge(){
		std::string name = e -> get_name();
		std::vector<std::string> names = { "HYDROGEN",\
						   "HELIUM",\
						   "LITHIUM",\
						   "BERYLLIUM",\
						   "BORON",\
						   "CARBON",\
						   "NITROGEN",\
						   "OXYGEN" };
		for ( int i = 0; i < names.size(); ++i ){
			if ( name == names[i] ){
				charge = i + 1;
			}
		}
	}

	void norm_atom(){
		for ( auto bf: e -> get_bf() ){
			bf -> renorm_bf();
		}
	}

	void new_coords( double x, double y, double z ){
		coords = _Coords(x,y,z);
	}

	int get_charge(){ return charge; }

private:
	int charge;
	_Coords coords;
	_Element * e;

};
