#pragma once

#include <iostream>  // стандартный ввод/вывод
#include <vector>    // библиотека для класса векторов

// Создадим класс для хранения трех натуральных чисел ---
// степеней полинома ( угловой части ):
// { i, j, k } --> x^i y^j z^k
class _Triple{

public:
	_Triple( int i, int j, int k ) {
		powers.push_back( i );
		powers.push_back( j );
		powers.push_back( k );
	}

	int get_i(){ return powers[0]; }
	int get_j(){ return powers[1]; }
	int get_k(){ return powers[2]; }

	void show(){
		std::cout << get_i() << " " \
			  << get_j() << " " \
			  << get_k() << std::endl;
	}

private:
	std::vector<int> powers;

};

// Введем функцию проверки равенства двух троек
bool triples_eq( _Triple t1, _Triple t2 ){
	if ( ( t1.get_i() != t2.get_i() ) || \
	   ( t1.get_j() != t2.get_j() ) || \
	   ( t1.get_k() != t2.get_k() ) ) { return false; } else { return true; }
}