#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "BasisFunctionClass.hpp"

// Создадим класс хим. элементов
// его объекты хранят имя элемента и вектор базисных функций
class _Element{

public:
	_Element( std::string name ): name(name) {} // конструктор
	~_Element() { // деструктор
		for( int i = 0; i < basis_functions.size(); ++i) {
			delete basis_functions[i];
		}
	}

// Добавим метод, дописывающий в вектор базисных функций новую функцию
	void add_basis_function( _Basis_function * bf ){
		basis_functions.push_back( bf );
	}

	int num_func(){
		int sum = 0;
		for ( auto bf: basis_functions ){
			sum += ( bf -> get_num_projections() );
		}
		return sum;
	}

	void show(){
		std::cout << "_______________________________________________________" << std::endl;
		std::cout << "!!! Successfully created new element: " << name << " !!!" << std::endl;
		std::cout << "_______________________________________________________" << std::endl;
		std::cout << std::endl;
	}

	void show_basis(){
		std::cout << "--> Basis consists of " << num_func() << " basis functions." << std::endl;
	}

	std::string get_name(){	return name; }

	std::vector<_Basis_function*> get_bf(){ return basis_functions; }

private:
	std::string name;
	std::vector<_Basis_function*> basis_functions; // вектор указателей на элементы класса _Basis_function
};

