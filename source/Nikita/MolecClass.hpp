#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <Eigen/Dense>

#include "CoordsClass.hpp"
#include "BasisClass.hpp"
#include "ElementClass.hpp"
#include "AtomClass.hpp"

class _Molecule{

public:
	_Molecule(){}
	~_Molecule(){}

	// Проверяем, что файл геометрии существует
	// если существует, вызываем функцию для его чтения
	void read_geom( std::string filename, _Basis * b ){

		std::ifstream fin( filename );

		if ( !fin ) throw std::invalid_argument( " Can't open a file " ); // если не существует, выдай ошибку
		else {
			std::cout << "File " << filename << " is opened" << std::endl; 
			parse_geom_file( fin, b ); // иначе --- читай его
		}

		fin.close();
	}


	// Читаем файл с геометрией
	void parse_geom_file( std::ifstream & fin, _Basis * b ){

		const int MAX_SIZE = 256;
		char line[ MAX_SIZE ];
		std::string current_string;
		_Atom * atom = NULL;

		for ( int i = 0; i < 2; ++i){
			fin.getline( line, MAX_SIZE );
		}

		while( fin.getline( line, MAX_SIZE ) ){
			current_string = line;

			std::stringstream ss( current_string );

			std::string name;
			double x, y, z;

			ss >> name >> x >> y >> z;

			for ( auto e: b -> get_elements() ){
				if ( name == ( e -> get_name() ) ){
					atom = new _Atom(x,y,z,e);
					atoms.push_back( atom );
				}
			}
		}
		show_geom();
		show_overlap_matrix();
	}

	void show_geom(){
		for ( auto a: atoms ){
			a -> show_atom();
		}
	}

	double E_coeff( int i, int j, int t, double Qx, double a, double b ){
		double p = a + b;
		double q = a * b / p;
		if ( ( t < 0 ) || ( t > (i + j) ) ) { return 0.0; }
		else if ( i == j == t == 0 ) { return exp( -q * Qx * Qx ); }
		else if ( j == 0 ){
			return ( 1.0 / (2.0 * p) ) * E_coeff (i-1, j, t-1, Qx, a, b) - \
			       ( q * Qx / a ) * E_coeff ( i-1, j, t, Qx, a, b) + \
			       ( t + 1 ) * E_coeff ( i-1, j, t+1, Qx, a, b );
		}
		else{
			return ( 1.0 / (2.0 * p) ) * E_coeff (i, j-1, t-1, Qx, a, b) + \
			       ( q * Qx / a ) * E_coeff ( i, j-1, t, Qx, a, b) + \
			       ( t + 1 ) * E_coeff ( i, j-1, t+1, Qx, a, b );
		}
	}

	double primitive_overlap( _Primitive p1, _Triple t1, _Coords c1, \
				  _Primitive p2, _Triple t2, _Coords c2 ){
		double S1 = E_coeff ( t1.get_i(), \
				      t2.get_i(), \
				      0, \
				      c1.get_x() - c2.get_x(), \
				      p1.get_alpha(), \
				      p2.get_alpha() );

		double S2 = E_coeff ( t1.get_j(), \
				      t2.get_j(), \
				      0, \
				      c1.get_y() - c2.get_y(), \
				      p1.get_alpha(), \
				      p2.get_alpha() );

		double S3 = E_coeff ( t1.get_k(), \
				      t2.get_k(), \
				      0, \
				      c1.get_z() - c2.get_z(), \
				      p1.get_alpha(), \
				      p2.get_alpha() );

		return S1 * S2 * S3 * pow( M_PI / ( p1.get_alpha() + p2.get_alpha() ), 1.5 );
	}

	double projection_overlap( _Projection p1, _Coords c1, _Projection p2, _Coords c2 ){
		double sum = 0.0;
		for ( auto i: p1.get_primitives() ){
			for ( auto j: p2.get_primitives() ){
				sum += i.get_coeff() * \
				       j.get_coeff() * \
				       primitive_overlap( i, p1.get_triple(), c1, j, p2.get_triple(), c2 );
			}
		}
		return sum;
	}

	
	Eigen::MatrixXd element_overlap( _Atom * a1, _Atom * a2 ){

		_Element * e1 = a1 -> get_e(), * e2 = a2 -> get_e();
		_Coords c1 = a1 -> get_c(), c2 = a2 -> get_c();

		Eigen::MatrixXd m( e1 -> num_func(), e2 -> num_func() );
		int i = 0;
		for ( auto bf1: e1 -> get_bf() ){
			for ( auto p1: bf1 -> get_projections() ){
				int j = 0;
				for ( auto bf2: e2 -> get_bf() ){
					for ( auto p2: bf2 -> get_projections() ){
						m(i, j) = projection_overlap( p1, c1, p2, c2 );
						++j;
					}
				}
				++i;
			}
		}

		return m;
	}

	int dimensions(){
		int sum = 0;
		for ( auto a: atoms ){
			sum += ( a -> get_e() ) ->  num_func();
		}
		return sum;
	}

	void show_overlap_matrix(){
		for ( int i = 0; i < atoms.size(); ++i ){
			for ( int j = 0; j < atoms.size(); ++j ){
				std::cout << "Overlap: " << atoms[i] -> get_e() -> get_name() << "_" << i+1 << " " \
							 << atoms[j] -> get_e() -> get_name() << "_" << j+1 <<  std::endl;
				std::cout << element_overlap( atoms[i], atoms[j] ) << std::endl;
				std::cout << std::endl;
			}
		}
	}

private:
	std::vector<_Atom*> atoms;
};
