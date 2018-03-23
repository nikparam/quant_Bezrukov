#pragma once

#include <vector>
#include <string>
#include "CGOBuilder.hpp"
#include "chemutils.hpp"

class Element
{
public:
    Element( std::string name ) : name(name)
    {
        charge = ChemUtils::GetZForAtom(name);
	}

	~Element()
	{
        for ( ContractedGaussianOrbital * cgo : CGOs )
            delete cgo;
	}

    void add_CGOs( CGOBuilder * CGObuilder );

	void show();

    ContractedGaussianOrbital * getCGO( const int n ) { return CGOs.at( n ); }

    size_t getCGOCount() { return CGOs.size(); }

    std::string getName() const { return name; }
    unsigned int getCharge() const { return charge; }

private:
    unsigned int charge;
	std::string name;
    std::vector<ContractedGaussianOrbital*> CGOs;
};

