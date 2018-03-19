#pragma once

#include <iostream>  // библиотека стандартного ввода/вывода
#include <cmath>     // библиотека математических функций
#include <vector>    // библиотека для класса векторов

#include "TripleClass.hpp"
#include "Misc.hpp"
#include "CoordsClass.hpp"

// Создадим класс, объекты которого --- тройки параметров гауссовых примитивов
// номер, множитель в экспоненте, коэффициент перед экспонентой
class _Primitive{

public:
	_Primitive( int num, double alpha, double coeff ): \
			num(num), alpha(alpha), coeff(coeff) { // инициализируем
	} 
	~_Primitive() { }

	// определим функции для вытаскивания свойств примитива (чтобы не делать их public переменными)
	int get_num( ) { return num; }
	double get_alpha( ) { return alpha; }
	double get_coeff( ) { return coeff; }

	// определим функцию для пересчета коэффициентов контрактации с учетом нормировки примитивов
	// в качестве входного параметра --- тройка показателей степеней полинома ( угловой части )
	void renorm( _Triple triple ){
		int i = triple.get_i();
		int j = triple.get_j();
		int k = triple.get_k();

		// Функция factor ( int i, double alpha ) лежит в Misc.hpp

		double N_x = factor( i, 2.0 * alpha );
		double N_y = factor( j, 2.0 * alpha );
		double N_z = factor( k, 2.0 * alpha );

		double N = std::sqrt( N_x * N_y * N_z );

		coeff /= N;
	}

private: // задаем внутренние переменные --- недоступны извне
	int num;
	double alpha, coeff;
};

