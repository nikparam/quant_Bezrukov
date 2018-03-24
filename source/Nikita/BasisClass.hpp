#pragma once

#include <iostream>  // библиотека стандартного ввода/вывода
#include <fstream>   // библиотека ввода/вывода из файла
#include <vector>    // библиотека для класса векторов
#include <string>    // библиотека для класса строк
#include <stdexcept> // библиотека для выдачи сообщений об ошибках
#include <locale>    // нужно в isdigit (?)
#include <sstream>   // для перевода данных из строки в переменные
#include <cmath>     // библиотека математических функций

#include "BasisFunctionClass.hpp"
#include "ElementClass.hpp"
#include "CoordsClass.hpp"

//Создадим класс базисных функций
//его объекты хранят базисные функции всех элементов, лежащих в файле filename
class _Basis{

public:
	_Basis( ) { } // конструктор
	~_Basis( ) {
		for( int i = 0; i < elements.size(); ++i ){ // деструктор
			delete elements[i];
		}
	}

	// Проверяем, что файл базиса существует
	// если существует, вызываем функцию для его чтения
	void read_basis( std::string filename ){

		std::ifstream fin( filename );

		if ( !fin ) throw std::invalid_argument( " Can't open a file " ); // если не существует, выдай ошибку
		else {
			std::cout << "File " << filename << " is opened" << std::endl; 
			parse_basis_file( fin ); // иначе --- читай его
		}

		fin.close();
	}

	//Функция для чтения файла
	// раскидывает информацию по соответствующим классам
	void parse_basis_file( std::ifstream & fin ){

		const int MAX_SIZE = 256;
		int primitives_num = 0;
		char line[ MAX_SIZE ];
		size_t first_ws, last_not_ws, num_primitives;
		std::string current_string, current_element;
		char first_char;
		_Element * element_pointer = NULL;
		_Basis_function * bf_pointer = NULL, * bf1_pointer = NULL;
		bool L_func = false;
		std::locale loc;

		// Будем читать файлик построчно

		while ( fin.getline( line, MAX_SIZE ) ){

			current_string = line;

		// Если строка начинается на !, $ или пустую строку --- не читaем ее, но проверяем, не пора ли добавить элемент в базис

			if ( current_string.size() == 0 || current_string.at(0) == '!' || current_string.at(0) == '$' ) {
				if ( element_pointer != NULL ) { // если указатель не пуст, то мы ранее проинициализировали элемент --- пора добавить его в базис
					elements.push_back( element_pointer ); // Добавляем в массив элементов новый элемент
//					show(); // выводим имя нового элемента
					element_pointer = NULL; // занулим, чтобы показать, что элемент добавлен
				} 
				continue;
			}

		// Ищем строчки, содержащие одно слово --- это названия наших элементов

			first_ws = current_string.find_first_of(' ');
			last_not_ws = current_string.find_last_not_of(' ');

		// Проверяем, что в строке первый пробел встречается только после всех слов --> мы нашли название элемента

			if ( ( first_ws == std::string::npos ) || ( first_ws > last_not_ws ) ) { 
				current_element = current_string.substr(current_string.find_first_not_of(' '), \
								        current_string.find_last_not_of(' ') + 1 ); // вытаскиваем имя
				element_pointer = new _Element( current_element ); // создаем указатель типа _Element
//				element_pointer -> show(); // выводим имя найденного элемента
				continue;
			}

		// Ищем строки с указанием симметрии угловых частей
		// Вытаскиваем первый значимый элемент строки

			first_char = current_string[ current_string.find_first_not_of(' ') ];

		// Если он не переводится в цифру, то это угловая часть

			if ( ! ( std::isdigit( first_char, loc ) ) ){
				if ( first_char != 'L' ){
					bf_pointer = new _Basis_function( first_char ); // создаем указатель типа _Basis_function
					primitives_num = std::stoi( &current_string[1] ); // определяем число примитивов
					continue;
				} else {
					L_func = true;
					bf_pointer = new _Basis_function( 'S' );
					bf1_pointer = new _Basis_function( 'P' );
					primitives_num = std::stoi( &current_string[1] );
					continue;
				}

			}else{ // иначе --- тройки параметров примитивов

				if ( !( L_func ) ){

					std::stringstream ss( current_string ); // задаем вывод строки в переменные
					int i;
					double alpha, coeff;

					ss >> i >> alpha >> coeff; // выводим значения в переменные
					bf_pointer -> add_primitive( i, alpha, coeff ); // по указателю добавляем новый примитив

					// если примитивы кончились, по указателю добавляем функцию

					if ( i == primitives_num ) {
						element_pointer -> add_basis_function( bf_pointer );
//						bf_pointer -> show_bf(); // выводим добавленную функцию
//						bf_pointer -> show_norm(); // выводим значение ее нормы ( если все хорошо, она равна 1 )
//						std::cout << std::endl;
					}
				} else {
					std::stringstream ss( current_string ); // задаем вывод строки в переменные
					int i;
					double alpha, coeff_s, coeff_p;

					ss >> i >> alpha >> coeff_s >> coeff_p; // выводим значения в переменные
					bf_pointer -> add_primitive( i, alpha, coeff_s ); // по указателю добавляем новый примитив
					bf1_pointer -> add_primitive( i, alpha, coeff_p ); // по указателю добавляем новый примитив

					if ( i == primitives_num ) {
						element_pointer -> add_basis_function( bf_pointer );
						element_pointer -> add_basis_function( bf1_pointer );
						L_func = false;
					}
				}
			}
		}
//		show_end(); // выводим сообщение о конце файла
	}

	std::vector<_Element*> get_elements() { return elements; }

	void show(){
		std::cout << std::endl;
		std::cout << "Complited basis with functions for: " << elements.end()[-1] -> get_name() << std::endl;
		std::cout << std::endl;
	}

	void show_end(){
		std::cout << "Basis file is scanned successfully" << std::endl;
		std::cout << "Proceed with the analisys" << std::endl;
		std::cout << std::endl;
	}

private:
	std::vector<_Element*> elements;
};


