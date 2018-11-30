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

    // для L функций
    CGOBuilder * CGObuilder_s = NULL;
    CGOBuilder * CGObuilder_p = NULL;

    bool isLfunction = false;

	// число примитивов в этой базисной функции
	int primitives_counter = 0; 

	while( infile.getline( buf, MAXLINE ) ) 
	{
		// текущая строчка в контейнере std::string
        current_string = buf;

        if ( DEBUG )
            std::cout << current_string << std::endl;

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

            if ( DEBUG ) std::cout << "Finalizing element." << std::endl;

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

            if ( DEBUG ) std::cout << "Creating new element: " << elp->getName() << std::endl;

            continue;
		}

		// если первый элемент в строчке это не число
		// то это строка вида "S  4"
		char first_symbol = current_string[current_string.find_first_not_of(' ')];

		if ( !std::isdigit( first_symbol, loc) )	
        {
            if ( first_symbol == 'L' )
            {
                isLfunction = true;

                // создаем S и P базисные функции
                CGObuilder_s = new CGOBuilder( 'S' );
                CGObuilder_p = new CGOBuilder( 'P' );

                if ( DEBUG ) std::cout << "Creating L function." << std::endl;
            }

            else
            {
                isLfunction = false;

                // создаем строителя базисных функций типа AngularPart (S, P, D и т.д.)
                CGObuilder = new CGOBuilder( first_symbol );

                // скармливаем функции stoi указатель на следующую за буквой позицию строки
                // например в строке "S  4" даем указатель на пробел, stoi убирает пробелы и возвращает (int) 4
                primitives_counter = std::stoi(&current_string[1]);

                if ( DEBUG ) std::cout << "Creating " << first_symbol << " function." << std::endl;
            }
        }

		// парсим строку вида "номер примитива, показатель экспоненты, коэффициент контрактации"
		// добавляем примитив в базисную функцию
		else
		{
            std::stringstream ss(current_string);

            int i; // номер примитива

            if ( isLfunction )
            {
                double alpha, coeff1, coeff2;
                // показатель экспоненты, коэффициент перед S функцией, коэффициент перед P функцией
                ss >> i >> alpha >> coeff1 >> coeff2;

                CGObuilder_s->add_primitive( alpha, coeff1 );
                CGObuilder_p->add_primitive( alpha, coeff2 );
            }
            else
            {
                double alpha, coeff;
                ss >> i >> alpha >> coeff;

                // билдер добавляет примитив ко всем базисным функциям, которые он строит
                CGObuilder->add_primitive( alpha, coeff );
            }

			// как только номер примитива совпадает с заявленным количеством примитивов, добавляем функцию текущему элементу
			if ( i == primitives_counter )
            {
                if ( isLfunction )
                {
                    //CGObuilder_s->normalize();
                    //CGObuilder_p->normalize();
                    elp->add_CGOs( CGObuilder_s );
                    elp->add_CGOs( CGObuilder_p );

                    if ( DEBUG ) std::cout << "Attaching L function to the element" << std::endl;

                    delete CGObuilder_s;
                    delete CGObuilder_p;
                }
                else
                {
                    //bfp->showQuantumNumbers();
                    // добавляем все построенные билдером базисные функции элементу
                    //CGObuilder->normalize();
                    elp->add_CGOs( CGObuilder );

                    if ( DEBUG ) std::cout << "Attaching CGO to the element." << std::endl;

                    // разрушаем строителя контрактированных орбиталей
                    delete CGObuilder;
                }
			}
		}
	}
}	

void Basis::show( const std::string type )
{
    if ( type == "short" )
    {
        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << "-------------- SHORT BASIS INFO ---------------" << std::endl;
        std::cout << "(Basis) No of elements: " << elements.size() << std::endl;
        for ( Element * el : elements )
            std::cout << "(Basis) Element: " << el->getName() << "; charge = " << el->getCharge() <<
                         "; number of CGOs: " << el->getCGOCount() << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;
    }

    if ( type == "full" )
    {
        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << "-------------- FULL BASIS INFO ----------------" << std::endl;
        std::cout << "No of elements: " << elements.size() << std::endl;
        for ( Element * el : elements )
        {
            el->show();
        }
        std::cout << "-----------------------------------------------" << std::endl;
    }
}
