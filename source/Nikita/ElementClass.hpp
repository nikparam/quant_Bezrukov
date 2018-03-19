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

	void add_coords( _Coords * c ){
		geom.push_back( c );
	}

	void show_geom(){
		std::cout << name << ": " \
			  << ( geom.end()[-1] -> get_x() ) << " " \
			  << ( geom.end()[-1] -> get_y() ) << " " \
			  << ( geom.end()[-1] -> get_z() ) << std::endl;
	}

	void show(){
		std::cout << "_______________________________________________________" << std::endl;
		std::cout << "!!! Successfully created new element: " << name << " !!!" << std::endl;
		std::cout << "_______________________________________________________" << std::endl;
		std::cout << std::endl;
	}

	std::string get_name(){	return name; }

	std::vector<_Basis_function*> get_bf(){ return basis_functions; }

private:
	std::string name;
	std::vector<_Basis_function*> basis_functions; // вектор указателей на элементы класса _Basis_function
	std::vector<_Coords*> geom;
};

