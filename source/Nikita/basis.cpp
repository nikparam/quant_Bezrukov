#include <iostream>  // библиотека стандартного ввода/вывода
#include <fstream>   // библиотека ввода/вывода из файла
#include <cmath>     // библиотека математических функций
#include <vector>    // библиотека для класса векторов
#include <string>    // библиотека для класса строк
#include <stdexcept> // библиотека для выдачи сообщений об ошибках
#include <locale>    // нужно в isdigit (?)
#include <sstream>   // для перевода данных из строки в переменные
#include <iomanip>   // для задания точности вывода

int d_factorial( int n ){
	return ( n <= 1 ) ? 1 : n * d_factorial( n - 2 );
}

double factor( int i, double alpha ){
	return  d_factorial( 2 * i - 1 ) * std::sqrt( M_PI ) / std::pow( 2, i ) / std::pow( alpha, i + 0.5 ) ;
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

bool triples_eq( _Triple t1, _Triple t2 ){
	if ( ( t1.get_i() != t2.get_i() ) || \
	   ( t1.get_j() != t2.get_j() ) || \
	   ( t1.get_k() != t2.get_k() ) ) { return false; } else { return true; }
}

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
	_Primitive( int num, double alpha, double coeff ): \
			num(num), alpha(alpha), coeff(coeff) { // инициализируем
	} 
	~_Primitive() { }

	// определим функции для вытаскивания свойств примитива (чтобы не делать их public переменными)
	int get_num( ) { return num; }
	double get_alpha( ) { return alpha; }
	double get_coeff( ) { return coeff; }

	// определим функцию для пересчета коэффициентов контрактации с учетом нормировки примитивов

	void renorm( _Triple triple ){
		int i = triple.get_i();
		int j = triple.get_j();
		int k = triple.get_k();

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

class _Projection{

public:
	_Projection( _Triple triple ): triple(triple) {}
	~_Projection(){}

	void add_primitive( int i, double alpha, double coeff ){
		primitives.emplace_back( i, alpha, coeff );
		primitives.end()[-1].renorm( triple );
	}

	double integral( _Primitive p1, _Primitive p2 ){

		int i = triple.get_i();
		int j = triple.get_j();
		int k = triple.get_k();

		double coeff1 = p1.get_coeff();
		double coeff2 = p2.get_coeff();

		double alpha1 = p1.get_alpha();
		double alpha2 = p2.get_alpha();

		double N_x = factor( i, alpha1 + alpha2 );
		double N_y = factor( j, alpha1 + alpha2 );
		double N_z = factor( k, alpha1 + alpha2 );

//		std::cout << "(" << i << j << k << ") " << \
//				    coeff1 * coeff2 * N_x * N_y * N_z << std::endl;
		return  coeff1 * coeff2 * N_x * N_y * N_z;
	}

	double count_norm(){
		double sum = 0.0;
		for ( int i = 0; i < primitives.size(); ++i ){
			for ( int j = 0; j < primitives.size(); ++j ){
				sum += integral ( primitives[i], primitives[j] ); 
			}
		}
		return sum; 
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

// Создадим класс базисных функций
// его элементы хранят символ угловой части и вектор примитивов
class _Basis_function{

public:
	_Basis_function( char angular_part ): angular_part( angular_part ) { add_triples(); } // инициализируем
	~_Basis_function(){ //деструктор
	}

	void add_projections( ){
		for ( auto & triple: triples){
			projections.emplace_back( triple );
		}
	}

// Добавим метод, дописывающий в набор примитивов новую гауссову функцию
// функция emplace_back записывает в конец вектора primitives объект класса примитив с параметрами i, alpha, coeff
	void add_primitive( int num, double alpha, double coeff ){
		for ( int i = 0; i < triples.size(); ++i ){
			projections[ i ].add_primitive( num, alpha, coeff );
		}
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

		add_projections();

	}

	char get_ap(){ return angular_part; }
	std::vector<_Triple> get_triples(){ return triples; }
	std::vector<_Projection> get_projections(){ return projections; }

	void show_bf(){
		std::cout << "--> Successfully created new basis function  of " << angular_part \
			  << " angular simmetry, contracted from: " << projections.end()[-1].get_primitives().size() \
			  << " primitive(s)." << std::endl;
		std::cout << std::endl;

	}

	void show_norm(){
		for ( auto & projection: projections ){
			std::cout << "(" << projection.get_triple().get_i() \
					 << projection.get_triple().get_j() \
					 << projection.get_triple().get_k() << "): N = "  << projection.get_primitives().size() << \
								               " norm = " << projection.count_norm() << std::endl;
			projection.show();
		}
	}

	void show_t(){
		for ( auto el : triples )
			el.show();
			std::cout << std::endl;
	}

private: // задаем параметры
	char angular_part;
	std::vector <_Triple> triples;
	std::vector <_Projection> projections;
};

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

// Добавим метод, дописывающий в вектор базовых функций новую функцию
	void add_basis_function( _Basis_function * bf ){
		basis_functions.push_back( bf );
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
					show();
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
				element_pointer -> show();
				continue;
			}

		// Ищем строки с указанием симметрии угловых частей
		// Вытаскиваем первый значимый элемент строки

			first_char = current_string[ current_string.find_first_not_of(' ') ];

		// Если он не переводится в цифру, то это угловая часть

			if ( ! ( std::isdigit( first_char, loc ) ) ){
				bf_pointer = new _Basis_function( first_char ); // создаем указатель типа _Basis_function
				primitives_num = std::stoi( &current_string[1] ); // определяем число примитивов
				continue;

			}else{ // иначе --- тройки параметров примитивов

				std::stringstream ss( current_string ); // задаем вывод строки в переменные
				int i;
				double alpha, coeff;

				ss >> i >> alpha >> coeff; // выводим значения в переменные
				bf_pointer -> add_primitive( i, alpha, coeff ); // по указателю добавляем новый примитив

				// если примитивы кончились, по указателю добавляем функцию

				if ( i == primitives_num ) {
					element_pointer -> add_basis_function( bf_pointer );
					bf_pointer -> show_bf();
					bf_pointer -> show_norm();
//					bf_pointer -> show_t();
					std::cout << std::endl;
				}
			}
		}
		show_end();
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
	std::cout << std::fixed << std::setprecision(7);

	std::string filename = "./tests/basis/cc-pvdz.gamess-us.dat";
	_Basis bs;
	bs.read( filename ); 
	
	return 0;

}
