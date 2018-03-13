#pragma once

#include <iostream>  // библиотека стандартного ввода/вывода
#include <cmath>     // библиотека математических функций
#include <vector>    // библиотека для класса векторов

// Создадим класс, объекты которого --- тройки параметров гауссовых примитивов
// номер, множитель в экспоненте, коэффициент перед экспонентой
class _Primitive{

public:
	_Primitive( int i, double alpha, double coeff, _Triple powers ): \
			i(i),     alpha(alpha), coeff(coeff),  powers(powers) { // инициализируем
		renormalize();
	} 
	~_Primitive() { }

	// определим функции для вытаскивания свойств примитива (чтобы не делать их public переменными)
	int get_i( ) { return i; }
	double get_alpha( ) { return alpha; }
	double get_coeff( ) { return coeff; }
	_Triple get_powers() { return powers; }

	// определим функцию для пересчета коэффициентов контрактации с учетом нормировки примитивов

	void renormalize( ){
		int i = powers.get_i();
		int j = powers.get_j();
		int k = powers.get_k();

		double N_x = factor( i, 2.0 * alpha );
		double N_y = factor( j, 2.0 * alpha );
		double N_z = factor( k, 2.0 * alpha );
		double N = std::sqrt( N_x * N_y * N_z );

		coeff /= N;
	}

private: // задаем внутренние переменные --- недоступны извне
	int i;
	double alpha, coeff;
	_Triple powers;
};


