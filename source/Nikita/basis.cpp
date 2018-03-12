#include <iostream>  // библиотека стандартного ввода/вывода
#include <fstream>   // библиотека ввода/вывода из файла
#include <cmath>     // библиотека математических функций
#include <vector>    // библиотека для класса векторов
#include <string>    // библиотека для класса строк
#include <stdexcept> // библиотека для выдачи сообщений об ошибках
#include <locale>    // нужно в isdigit (?)
#include <sstream>   // для перевода данных из строки в переменные
#include <iomanip>   // для задания точности вывода

int factorial( int n ){

	return ( n == 0 || n == 1) ? 1 : n * factorial( n - 1 );
}

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

class _Center{

public:
	_Center( double x, double y, double z ){
		coords.push_back( x );
		coords.push_back( y );
		coords.push_back( z );
	}
private:
	std::vector<double> coords;

};


// Создадим класс, объекты которого --- тройки параметров гауссовых примитивов
// номер, множитель в экспоненте, коэффициент перед экспонентой
class _Primitive{

public:
	_Primitive( int i, double alpha, double coeff ): i(i), alpha(alpha), coeff(coeff) {} // инициализируем
	~_Primitive() {}

	// определим функции для вытаскивания свойств примитива (чтобы не делать их public переменными)
	int get_i( ) { return i; }
	double get_alpha( ) { return alpha; }
	double get_coeff( ) { return coeff; }

	// определим функцию для пересчета коэффициентов контрактации с учетом нормировки примитивов
	double renorm_factor( int n ){
		return factorial( n ) / std::pow( 2, 2 * n ) * std::sqrt( M_PI / ( 2.0 * alpha ) ) /  std::pow( alpha, n );
	}

private: // задаем внутренние переменные --- недоступны извне
	int i;
	double alpha, coeff;
};


// Создадим класс базисных функций
// его элементы хранят символ угловой части и вектор примитивов
class _Basis_function{

public:
	_Basis_function( char angular_part ): angular_part( angular_part ) {} // инициализируем
	~_Basis_function(){ //деструктор
		for( int i = 0; i < primitives.size(); ++i) primitives[i].~_Primitive();
	}


// Добавим метод, дописывающий в набор примитивов новую гауссову функцию
// функция emplace_back записывает в конец вектора primitives объект класса примитив с параметрами i, alpha, coeff
	void add_primitive( int i, double alpha, double coeff ){
		primitives.emplace_back( i, alpha, coeff );
	}

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
					if ( i + j + k == l ) triples.emplace_back( i, j, k );
				}
			}
		}

	}

/*	void renormalization(){

		for( auto a : primitives ){

			for( auto t: triples ){

				double N_x = a.renorm_factor( t.get_i() );
				double N_y = a.renorm_factor( t.get_j() );
				double N_z = a.renorm_factor( t.get_k() );

			 	std::cout << angular_part << " " << std::sqrt( N_x * N_y * N_z ) << std::endl;

			}
		}
	} */

	char get_ap(){ return angular_part; }

	std::vector<_Primitive> get_primitives(){ return primitives; }

	std::vector<_Triple> get_triples(){ return triples; }

	void show_bf(){
		std::cout << "--> Successfully created new basis function  of " << angular_part \
			  << " angular simmetry, contracted from: " << primitives.size() \
			  << " primitive(s)." << std::endl;
		std::cout << std::endl;

	}
	void show_p(){
		std::cout << primitives.end()[-1].get_i() << ") Primitive "  \
			  << " with exp. mult. " << primitives.end()[-1].get_alpha() \
			  << " and contract. coeff. " << primitives.end()[-1].get_coeff() << std::endl;
	}
	void show_t(){
		for ( auto el : triples )
			el.show();
			std::cout << std::endl;
	}

private: // задаем параметры
	char angular_part;
	std::vector <_Primitive> primitives;
	std::vector <_Triple> triples;
};

// Создадим класс хим. элементов
// его объекты хранят имя элемента и вектор базисных функций
class _Element{

public:
	_Element( std::string name ): name(name) {} // конструктор
	~_Element() { // деструктор
		for( int i = 0; i < basis_functions.size(); ++i) {
			delete basis_functions[i];
//			basis_functions[i] -> ~_Basis_function();
		}
	}

// Добавим метод, дописывающий в вектор базовых функций новую функцию
	void add_basis_function( _Basis_function * bf ){
		basis_functions.push_back( bf );
	}

	void norm(){

		double sum = 0.0;

		for ( auto bf : basis_functions ){
			for ( auto t: bf -> get_triples() ){
				for ( auto p: bf -> get_primitives() ){
					double N_x = p.renorm_factor( t.get_i() );
					double N_y = p.renorm_factor( t.get_j() );
					double N_z = p.renorm_factor( t.get_k() );

					sum += p.get_coeff() * p.get_coeff() * N_x * N_y * N_z ;
				}
			}
		}

		std::cout << name << ": " << sum << std::endl;

	}

	void show(){
		std::cout << "_______________________________________________________" << std::endl;
		std::cout << "!!! Successfully created new element: " << name << " !!!" << std::endl;
		std::cout << "_______________________________________________________" << std::endl;
		std::cout << std::endl;
	}

	std::string get_name(){
		return name;
	}

private:// инициализируем
	std::string name;
	std::vector<_Basis_function*> basis_functions; // вектор указателей на элементы класса _Basis_function

};

//Создадим класс базисных функций
//его объекты хранят базисные функции всех элементов, лежащих в файле filename
class _Basis{

public:
	_Basis( ) { } // конструктор
	~_Basis( ) {
		for( int i = 0; i < elements.size(); ++i ){ // деструктор
			delete elements[i];
//			elements[i] -> ~_Element();
		}
	}

// Проверяем, что файл существует
// если существуем, вызываем функцию для его чтения
	void read( std::string filename ){

		std::ifstream fin( filename );

		if ( !fin ) throw std::invalid_argument( " Can't open a file " ); // если не существует, выдай ошибку
		else {
			std::cout << "File is opened" << std::endl; 
			parse_file( fin ); // иначе --- читай его
		}

		fin.close();
	}

//Функция для чтения файла
// раскидывает информацию по соответствующим классам
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

		// Будем читать файлик построчно

		while ( fin.getline( line, MAX_SIZE ) ){

			current_string = line;

		// Если строка начинается на !, $ или пустую строку --- не читaем ее, но проверяем, не пора ли добавить элемент в базис

			if ( current_string.size() == 0 || current_string.at(0) == '!' || current_string.at(0) == '$' ) {
				if ( element_pointer != NULL ) { // если указатель не пуст, то мы ранее проинициализировали элемент --- пора добавить его в базис
					elements.push_back( element_pointer ); // Добавляем в массив элементов новый элемент
					element_pointer -> norm();
//					show();
					element_pointer = NULL; // занулим, чтобы показать, что элемент добавлен
				} 
				continue;
			}

		// Ищем строчки, содержащие одно слово --- это названия наших элементов

			first_ws = current_string.find_first_of(' ');
			last_not_ws = current_string.find_last_not_of(' ');

		// Проверяем, что в строке первый пробел встречается только после всех слов

			if ( ( first_ws == std::string::npos ) || ( first_ws > last_not_ws ) ) { 
				current_element = current_string.substr(current_string.find_first_not_of(' '), \
								        current_string.find_last_not_of(' ') + 1 ); // вытаскиваем имя
				element_pointer = new _Element( current_element ); // создаем указатель типа _Element
//				element_pointer -> show();
				continue;
			}

		// Ищем строки с указанием симметрии угловых частей
		// Вытаскиваем первый значимый элемент строки

			first_char = current_string[ current_string.find_first_not_of(' ') ];

		// Если он не переводится в цифру, то это угловая часть

			if ( ! ( std::isdigit( first_char, loc ) ) ){
				bf_pointer = new _Basis_function( first_char ); // создаем указатель типа _Basis_function
				bf_pointer -> add_triples();
//				bf_pointer -> show_t();
				primitives_num = std::stoi( &current_string[1] ); // определяем число примитивов
				continue;

			}else{ // иначе --- тройки параметров примитивов

				std::stringstream ss( current_string ); // задаем вывод строки в переменные
				int i;
				double alpha, coeff;

				ss >> i >> alpha >> coeff; // выводим значения в переменные
				bf_pointer -> add_primitive( i, alpha, coeff ); // по указателю добавляем новый примитив
//				bf_pointer -> show_p(); 

				// если примитивы кончились, по указателю добавляем функцию

				if ( i == primitives_num ) {
					element_pointer -> add_basis_function( bf_pointer );
//					bf_pointer -> renormalization();
//					bf_pointer -> show_bf(); 
				}
			}
		}
//		show_end();
	}

	void show(){
		std::cout << std::endl;
		std::cout << "Complited basis with functions for: " << elements.end()[-1] -> get_name() << std::endl;
		std::cout << std::endl;
	}

	void show_end(){
		std::cout << "Basis file is scanned successfully" << std::endl;
		std::cout << "Proceed with the analisys" << std::endl;
	}

private:
	std::vector<_Element*> elements;
};


int main(){
	std::cout << std::fixed << std::setprecision(6);

	std::string filename = "./tests/basis/cc-pvdz.gamess-us.dat";
	_Basis bs;
	bs.read( filename ); 
	
	return 0;

}
