#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <locale>
#include <sstream>
#include <iomanip>

// Создадим класс, объекты которого --- тройки параметров гауссовых примитивов
// номер, множитель в экспоненте, коэффициент перед экспонентой
class _Primitive{

public:
	_Primitive( int i, double alpha, double coeff ): i(i), alpha(alpha), coeff(coeff) {} // инициализируем

	int get_i( ) { return i; }
	double get_alpha( ) { return alpha; }
	double get_coeff( ) { return coeff; }		
private: // задаем внутренние переменные --- недоступны изве
	int i;
	double alpha, coeff;
};


// Создадим класс базисных функций
// его элементы хранят символ угловой части и вектор примитивов
class _Basis_function{

public:
	_Basis_function( char angular_part ): angular_part( angular_part ) {} // инициализируем

// Добавим метод:
// функция emplace_back записывает в конец вектора primitives объект класса примитив с параметрами i, alpha, coeff
	void add_primitive( int i, double alpha, double coeff ){
		primitives.emplace_back( i, alpha, coeff );
	}

	void show_bf(){
		std::cout << "Successfully created new basis function  of " << angular_part \
			  << " angular simmetry, contracted from: " << primitives.size() \
			  << " primitives." << std::endl;

	}
	void show_p(){
		std::cout << "Primitive: " <<  primitives.end()[-1].get_i() \
			  << " " << primitives.end()[-1].get_alpha() \
			  << " " << primitives.end()[-1].get_coeff() << std::endl;
	}

private: // задаем параметры
	char angular_part;
	std::vector <_Primitive> primitives;
};

// Создадим класс хим. элементов
// его объекты хранят имя элемента и вектор базисных функций
class _Element{

public:
	_Element( std::string name ): name(name) {} // инициализируем

// Добавим метод
	void add_basis_function( _Basis_function & bf ){
		basis_functions.emplace_back( bf );
	}

	void show(){
		std::cout << "Successfully created new element: " << name << std::endl;
	}

private:
	std::string name;
	std::vector<_Basis_function*> basis_functions;

};

class _Basis{

public:
	_Basis( ) { }
	~_Basis( ) { }	

	void read( std::string filename ){

		std::ifstream fin( filename );

		if ( !fin ) throw std::invalid_argument( " Can't open a file " );
		else {
			std::cout << "File is opened" << std::endl;
			parse_file( fin );
		}

		fin.close();
	}

	void parse_file( std::ifstream & fin ){

		const int MAX_SIZE = 256;
		int primitives_num = 0;
		char line[ MAX_SIZE ];
		size_t first_ws, last_not_ws, num_primitives;
		std::string current_string, current_element;
		char first_char;
		_Element * element_pointer = NULL;
		_Basis_function * bf_pointer = NULL;
		std::locale loc;		

		// Будем парсить файлик
		while ( fin.getline( line, MAX_SIZE ) ){
			current_string = line;
		// Если строка начинается на ! или $ --- не чиатем ее
			if ( current_string.size() == 0 || current_string.at(0) == '!' || current_string.at(0) == '$' ) continue;

		// В противном случае 
		// Ищем строчки, начинающиеся со слова --- это наши элементы
			first_ws = current_string.find_first_of(' ');
			last_not_ws = current_string.find_last_not_of(' ');

		// Проверяем, что в строке первый пробел встречается после всех слов
			if ( ( first_ws == std::string::npos ) || ( first_ws > last_not_ws ) ) { 
				current_element = current_string.substr(current_string.find_first_not_of(' '), \
								        current_string.find_last_not_of(' ') + 1 ); // вытаскиваем имя
				element_pointer = new _Element( current_element ); // создаем указатель типа _Element
				element_pointer -> show();
				continue;
			}

		// Вытаскиваем первый значимый элемент строки
			first_char = current_string[ current_string.find_first_not_of(' ') ];
		// Если он не переводится в цифру, то это угловая часть
			if ( ! ( std::isdigit( first_char, loc ) ) ){
				bf_pointer = new _Basis_function( first_char ); // создаем указатель типа _Basis_function
				primitives_num = std::stoi( &current_string[1] ); // определяем число примитивов
				bf_pointer -> show_bf();
				continue;
			}else{
				std::stringstream ss( current_string );
				int i;
				double alpha, coeff;
				ss >> i >> alpha >> coeff;
				bf_pointer -> add_primitive( i, alpha, coeff );
				bf_pointer -> show_p();
				if ( i == primitives_num ) { element_pointer -> add_basis_function( bf_pointer ) }
			}
		}
	}

private:

};


int main(){
	std::cout << std::fixed << std::setprecision(6);

	std::string filename = "./tests/basis/cc-pvdz.gamess-us.dat";
	_Basis bs;
	bs.read( filename ); 
	return 0;

}
