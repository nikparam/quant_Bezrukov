#pragma once

#include <vector>
#include <string>

using std::vector;
using std::string;

class Element
{
public:
	Element( string name ) : name(name)
	{
	}

	~Element()
	{
		for ( BasisFunction * bf : basis_functions )
			delete bf;
	}

	void add_basis_function( BasisFunction * bf )
	{
		basis_functions.push_back( bf );
	}

	void show()
	{
		cout << "Element name: " << name << endl;
	}

private:
	string name;
	vector<BasisFunction*> basis_functions;
};

