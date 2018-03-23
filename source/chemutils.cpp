#include "chemutils.hpp"

namespace ChemUtils
{
    const char * atoms[] =
    {
        "HYDROGEN",
        "HELIUM",
        "LITHIUM",
        "BERYLLIUM",
        "BORON",
        "CARBON",
        "NITROGEN",
        "OXYGEN",
        "FLUORINE",
        "NEON"
    };


    unsigned int GetZForAtom( std::string & name )
    {
        int size = 10;
        for ( int i = 0; i < size; i++ )
            if ( name == atoms[i])
                return i + 1;

        std::invalid_argument("Unknown element!");
    }
}
