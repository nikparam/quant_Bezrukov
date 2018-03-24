#include <iostream>  // библиотека стандартного ввода/вывода
#include <fstream>   // библиотека ввода/вывода из файла
#include <cmath>     // библиотека математических функций
#include <vector>    // библиотека для класса векторов
#include <string>    // библиотека для класса строк
#include <stdexcept> // библиотека для выдачи сообщений об ошибках
#include <locale>    // нужно в isdigit (?)
#include <sstream>   // для перевода данных из строки в переменные
#include <iomanip>   // для задания точности вывода

#include "Misc.hpp"
#include "TripleClass.hpp"
#include "PrimitiveClass.hpp"
#include "ProjectionClass.hpp"
#include "BasisFunctionClass.hpp"
#include "ElementClass.hpp"
#include "BasisClass.hpp"
#include "CoordsClass.hpp"
#include "MolecClass.hpp"

int main(){
	std::cout << std::fixed << std::setprecision(7);

	std::string basis_name;
	std::cin >> basis_name;
	_Basis bs;
	bs.read_basis( basis_name );

	std::string geom_name;
	std::cin >> geom_name;
	_Molecule m;
	m.read_geom( geom_name, & bs );
	
	return 0;

}
