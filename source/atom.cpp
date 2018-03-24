#include "atom.h"

Atom::Atom( std::string name, double x, double y, double z)
    : name(name), x(x), y(y), z(z)
{
}

void Atom::attachElementFunctions( Basis * basis )
{
    for ( size_t i = 0; i < basis->getElementsCount(); i++ )
    {
        if ( name == basis->getElement(i)->getName() )
            element = basis->getElement(i);
    }
}
