#pragma once

#include <iostream>
#include <cmath>
#include <vector>

#include "TripleClass.hpp"
#include "PrimitiveClass.hpp"
#include "Misc.hpp"

// Создадим класс, в котором будет хранится одна базисная функция --- вектор примитивов, 
// отвечающий одной и той же тройке показателей полинома ( угловой части )
class _Projection{

public:
	_Projection( _Triple triple ): triple(triple) {}
	~_Projection(){}

	// Создадим функцию для заполнения массива примитивов
	// в ходе заполнения новый примитив сразу нормируется
	void add_primitive( int i, double alpha, double coeff ){
		primitives.emplace_back( i, alpha, coeff );
		primitives.end()[-1].renorm( triple );
	}

	// Создадим функцию для расчета интеграла двух примитивов
	double integral( _Primitive p1, _Primitive p2 ){

		int i = triple.get_i();
		int j = triple.get_j();
		int k = triple.get_k();

		double coeff1 = p1.get_coeff();
		double coeff2 = p2.get_coeff();

		double alpha1 = p1.get_alpha();
		double alpha2 = p2.get_alpha();

		// функция factor ( int i, double alpha ) лежит в Misc.hpp
		double N_x = factor( i, alpha1 + alpha2 );
		double N_y = factor( j, alpha1 + alpha2 );
		double N_z = factor( k, alpha1 + alpha2 );

		return  coeff1 * coeff2 * N_x * N_y * N_z;
	}

	// Создадим функцию для расчета нормы базисной функции
	double count_norm(){
		double sum = 0.0;
		for ( int i = 0; i < primitives.size(); ++i ){
			for ( int j = 0; j < primitives.size(); ++j ){
				sum += integral ( primitives[i], primitives[j] ); 
			}
		}
		return sum; 
	}

	void renorm_projection(){
		double N = count_norm();
		for ( auto & p : primitives ){
			p.total_renorm( N );
		}
	}

	void show(){
		for( auto & primitive: primitives ){
			std::cout << primitive.get_num() << ") alpha = " << primitive.get_alpha() \
					               << " coeff = "  << primitive.get_coeff() << std::endl;
		}
	}

	_Triple get_triple() { return triple; }
	std::vector<_Primitive> get_primitives() { return primitives; }

private:
	_Triple triple;
	std::vector<_Primitive> primitives;

};

