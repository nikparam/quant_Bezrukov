#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "Basis.hpp"
#include "Element.hpp"

class Atom
{
public:
    Atom( std::string name, double x, double y, double z );

    void attachElementFunctions( Basis * basis );

    double get_x() const { return x; }
    double get_y() const { return y; }
    double get_z() const { return z; }
    std::string get_name() const { return name; }

    std::vector<double> getCoords() const {
        std::vector<double> temp{x, y, z};
        return temp;
    }

    Element * get_element() const { return element; }

private:
    std::string name;
    Element * element;

    double x;
    double y;
    double z;
};
