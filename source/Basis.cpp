#include "Basis.hpp"

void Basis::read( std::string filename )
{
    std::ifstream infile( filename );
	
    if ( !infile )
		throw std::invalid_argument("Can't open the file");
	else
		parse_file( infile );

	infile.close();	
}

void Basis::parse_file( std::ifstream & infile )
{
	std::locale loc;

	std::string current_string;
	std::string current_element;

	char buf[ MAXLINE ];
	
	// указатель на элемент, который мы заполняем
	Element * elp = NULL;
    // указатель на текущий строитель CGO (контрактированных гауссовых орбиталей)
    CGOBuilder * CGObuilder = NULL;

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
		if ( first_ws == std::string::npos || (first_ws > last_not_ws) )
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
            // создаем строителя базисных функций типа AngularPart (S, P, D и т.д.)
            CGObuilder = new CGOBuilder( first_symbol );

             // скармливаем функции stoi указатель на следующую за буквой позицию строки
             // например в строке "S  4" даем указатель на пробел, stoi убирает пробелы и возвращает (int) 4
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

            // билдер добавляет примитив ко всем базисным функциям, которые он строит
            CGObuilder->add_primitive( alpha, coeff );

			// как только номер примитива совпадает с заявленным количеством примитивов, добавляем функцию текущему элементу
			if ( i == primitives_counter )
			{
                //bfp->showQuantumNumbers();
				// добавляем все построенные билдером базисные функции элементу
                elp->add_CGOs( CGObuilder );

                // разрушаем строителя контрактированных орбиталей
                delete CGObuilder;
			}
		}
	}
}	

void Basis::show( const std::string type )
{
    if ( type == "short" )
    {
        std::cout << "No of elements: " << elements.size() << std::endl;
        for ( Element * el : elements )
            std::cout << "Element: " << el->getName() << "; charge = " << el->getCharge() <<
                         "; number of CGOs: " << el->getCGOCount() << std::endl;
    }

    if ( type == "full" )
    {
        std::cout << "No of elements: " << elements.size() << std::endl;
        for ( Element * el : elements )
        {
            el->show();
        }
    }
}
