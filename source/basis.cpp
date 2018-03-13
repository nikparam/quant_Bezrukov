#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream> 
#include <sstream> // streamstring
#include <cmath>

#include "Vector3D.hpp"
#include "quantumNumbers.hpp"
#include "Primitive.hpp"
#include "BasisFunction.hpp"
#include "Element.hpp"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::invalid_argument;
using std::vector;

#define MAXLINE 100

class Basis
{
public:
	Basis( ) { }
	~Basis()
	{
		for ( Element * el : elements )
			delete el;
	}

	void read( string filename )
	{
		ifstream infile( filename );
		
		if ( !infile )
			throw invalid_argument("Can't open the file");
		else
			parse_file( infile );

		infile.close();	
	}
		
	void parse_file( ifstream & infile )
	{
		std::locale loc;

		string current_string;
    	string current_element;

		char buf[ MAXLINE ];
	
		string current_group = "";

		// указатель на элемент, который мы заполняем
		Element * elp = NULL;
		// указатель на текущий базисную функцию
		BasisFunction * bfp = NULL;
		
		// число примитивов в этой базисной функции
		int primitives_counter = 0; 

		while( infile.getline( buf, MAXLINE ) ) 
		{
			// текущая строчка в контейнере std::string
			current_string = buf;
	
			// пустая строчка или строчка с долларом ($END) являются признаками окончания базисного набора
			// текущего элемента. поэтому если мы встречаем такую строчку, если мы заполняли до этого базисный 
			// набор некоторого элемента, то мы запихиваем указатель на текущий элемент в наш контейнер элементов, 
			// а текущий указатель обнуляем.	
			if ( current_string.size() == 0 || current_string.at(0) == '!' || current_string.at(0) == '$' )
			{
				if ( elp != NULL )
				{
					elements.push_back( elp );
					elp = NULL;
				}

				continue;
			}	

			// позиция первого пробела
			size_t first_ws = current_string.find_first_of(" ");
			// позиция последнего не пробела
			size_t last_not_ws = current_string.find_last_not_of(" ");

			// если нет пробелов или пробелы идут в конце строки,
			// то значит у нас только одно слово, значит это название элемента
   			if ( first_ws == string::npos || (first_ws > last_not_ws) )
			{
				// обдираем строку от пробелов 
				current_element = current_string.substr(current_string.find_first_not_of(" "), current_string.find_last_not_of(" ") + 1);
				// создали новый элемент
				elp = new Element( current_element );

				continue;
			}

			// если первый элемент в строчке это не число
			// то это строка вида "S  4"
			char first_symbol = current_string[current_string.find_first_not_of(' ')];
			
			if ( !std::isdigit( first_symbol, loc) )
			{
				// создается базисная функция;
				// внутри BasisFunction генерируется набор троек QuantumNumbers, с которыми будут генерироваться примитивы внутри этой базисной функции
				bfp = new BasisFunction(first_symbol);
				
				// скармливаем функции stoi указатель на следующую за буквой позицию строки
				// например в строке "S  4" даем указатель на пробел, stoi стрипит пробелы и возвращает (int) 4
				primitives_counter = std::stoi(&current_string[1]);
			}

			// парсим строку вида "номер примитива, показатель экспоненты, коэффициент контрактации"
			// добавляем примитив в базисную функцию
			else
			{
				std::stringstream ss(current_string);
				int i;
				double alpha, coeff;
				ss >> i >> alpha >> coeff;
				bfp->add_primitive( alpha, coeff );

				// как только номер примитива совпадает с заявленным количеством примитивов, добавляем функцию текущему элементу
				if ( i == primitives_counter )
				{
					bfp->fillQuantumNumbers();
					//bfp->showQuantumNumbers();
					elp->add_basis_function( bfp ); 
				}
			}

		}
	}	

	void show()
	{
		cout << "no of elements: " << elements.size() << endl;
		for ( Element * el : elements )
		{
			el->show();	
		}
	}

private:
	vector<Element*> elements;
};


int main()
{
	Basis basis;
	basis.read("../tests/basis/cc-pvdz.gamess-us.dat");

	basis.show();

	return 0;
}

