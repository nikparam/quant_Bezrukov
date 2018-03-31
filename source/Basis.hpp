#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept> // invalid_argument
#include <fstream> 
#include <sstream> // streamstring
#include <cmath>

#include "QuantumNumbers.hpp"
#include "Primitive.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Element.hpp"

class Basis
{
public:
	Basis( ) { }
	~Basis()
	{
		for ( Element * el : elements )
			delete el;
	}

	void read( std::string filename );
	void parse_file( std::ifstream & infile );

    void show( const std::string type );

	// at -- проверяет, что элемент с таким номером находится в векторе
	Element * getElement( const int n ) { return elements.at(n); }
	size_t getElementsCount() { return elements.size(); }

private:
    const int MAXLINE = 200;
    vector<Element*> elements;
};

