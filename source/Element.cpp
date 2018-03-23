#include "Element.hpp"

void Element::add_CGOs( CGOBuilder * CGObuilder )
{
    for ( ContractedGaussianOrbital * CGO : CGObuilder->getCGOs() )
        CGOs.push_back( CGO );
}


void Element::show()
{
    std::cout << "Element name: " << name << "; charge = " << charge << std::endl;
    for ( ContractedGaussianOrbital * CGO : CGOs )
        CGO->show();
}

