#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "TripleClass.hpp"
#include "ProjectionClass.hpp"

// Создадим класс базисных функций
// его элементы хранят символ угловой части --> вектор троек степеней полиномов ( угловых частей )
// и вектор проекций ( их столько, сколько позволяет максимальный орбитальный момент )
class _Basis_function{

public:
	_Basis_function( char angular_part ): angular_part( angular_part ) { add_triples(); } // инициализируем
	~_Basis_function(){ //деструктор
	}

	// Создадим функцию для заполнения массива проекций
	// их создается столько же, сколько возможно троек { i, j, k }
	void add_projections( ){
 		for ( auto triple: triples){
			_Projection * p = new _Projection( * triple );
			projections.push_back( p );
		}
	}

	// Создадим функцию, которая к каждой проекции добавляет новый примитив
	void add_primitive( int num, double alpha, double coeff ){
		for ( int i = 0; i < triples.size(); ++i ){
			projections[ i ] -> add_primitive( num, alpha, coeff );
		}
	}

	// Создадим функцию, которая по символу угловой части генерирует массив троек { i, j, k }
	void add_triples(){

		int l;

		switch ( angular_part ){

			case 'S':
			{
				l = 0;
				break;
			}

			case 'P':
			{
				l = 1;
				break;
			}

			case 'D':
			{
				l = 2;
				break;
			}

			case 'F':
			{
				l = 3;
				break;
			}

			default:
			{
				throw std::invalid_argument("Can't identify angular part");
			}
		}

		for( int i = 0; i <= l; ++i ){
			for( int j = 0; j <= l - i; ++j ){
				for( int k = 0; k <= l - i - j; ++k ){
					if ( i + j + k == l ) { 
						_Triple * t = new _Triple( k, j, i );
						triples.push_back( t );
					}
				}
			}
		}

		add_projections();

	}

	void renorm_bf(){
		for ( auto pr: projections ){
			std::cout << pr -> count_norm();
			pr -> renorm_projection();
			std::cout << " " << pr -> count_norm() << std::endl;
		}
	}

	char get_ap(){ return angular_part; }
	std::vector<_Triple*> get_triples(){ return triples; }
	std::vector<_Projection*> get_projections(){ return projections; }
	int get_num_projections(){ return projections.size(); }

	void show_bf(){
		std::cout << "--> Successfully created new basis function  of " << angular_part \
			  << " angular simmetry, contracted from: " << projections.end()[-1] -> get_primitives().size() \
			  << " primitive(s)." << std::endl;
		std::cout << std::endl;

	}

	void show_norm(){
		for ( auto projection: projections ){
			std::cout << "(" << projection -> get_triple().get_i() \
					 << projection -> get_triple().get_j() \
					 << projection -> get_triple().get_k() << \
				     "): N = "  << projection -> get_primitives().size() << \
				     " norm = " << projection -> count_norm() << std::endl;
			projection -> show();
		}
	}

	void show_t(){
		for ( auto el : triples )
			el -> show();
			std::cout << std::endl;
	}

private: // задаем параметры
	char angular_part;
	std::vector <_Triple*> triples;
	std::vector <_Projection*> projections;
};

