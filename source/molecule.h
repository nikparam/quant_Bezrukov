#pragma once

#include "atom.h"
#include "Basis.hpp"
#include "ContractedGaussianOrbital.hpp"
#include <iomanip>

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>

class Molecule
{
public:
    Molecule();

    void readGeometryFile( std::string const & filename );
    void parseGeometryFile( std::ifstream & infile );

    double calculateEijt( int i, int j, int t, double Qx, double a, double b );
    double overlapPrimitive( Primitive * a, Primitive * b,
                             QuantumNumbers * t1, QuantumNumbers * t2,
                             std::vector<double> coords1,
                             std::vector<double> coords2 );

    double overlapCGO( ContractedGaussianOrbital * cgo1, ContractedGaussianOrbital * cgo2,
                       std::vector<double> coords1, std::vector<double>  coords2 );

    void fillOverlapMatrix( );
    void saveOverlapMatrix( std::string filename = "" );

    Atom * getAtom( int n ) { return &atoms.at(n); }

    void setBasis( Basis * basis ) { this->basis = basis; }
    void showAtoms();

private:
    const int MAXLINE = 100;

    std::vector<Atom> atoms;
    Basis * basis;

    Eigen::MatrixXd overlapMatrix;
};
