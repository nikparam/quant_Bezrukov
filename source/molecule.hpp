#pragma once

#include "atom.hpp"
#include "Basis.hpp"
#include "ContractedGaussianOrbital.hpp"
#include <iomanip>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

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

    std::vector<double> gaussianProductCenter( double a, std::vector<double> & A,
                                               double b, std::vector<double> & B );

    double calculateEijt( int i, int j, int t, double Qx, double a, double b );
    double calculateHCintegral( int t, int u, int v, int n, double p, double PCx, double PCy, double PCz, double RPC );

    double nuclearAttractionPrimitive( Primitive * a, QuantumNumbers * t1, std::vector<double> coords1,
                                       Primitive * b, QuantumNumbers * t2, std::vector<double> coords2,
                                       std::vector<double> nuclearCoords );

    double nuclearAttractionCGO( ContractedGaussianOrbital * cgo1, ContractedGaussianOrbital * cgo2,
                                 std::vector<double> coords1, std::vector<double> coords2,
                                 std::vector<double> nuclearCoords );

    double nuclearAttractionCGO_total( ContractedGaussianOrbital * cgo1, ContractedGaussianOrbital * cgo2,
                                       std::vector<double> coords1, std::vector<double> coords2 );

    double overlapPrimitive( Primitive * a, Primitive * b,
                             QuantumNumbers * t1, QuantumNumbers * t2,
                             std::vector<double> coords1,
                             std::vector<double> coords2 );

    double overlapCGO( ContractedGaussianOrbital * cgo1, ContractedGaussianOrbital * cgo2,
                       std::vector<double> coords1, std::vector<double>  coords2 );

    double kineticPrimitive( Primitive * a, Primitive * b,
                             QuantumNumbers * t1, QuantumNumbers * t2,
                             std::vector<double> coords1,
                             std::vector<double> coords2 );
    double kineticCGO( ContractedGaussianOrbital * cgo1, ContractedGaussianOrbital * cgo2,
                       std::vector<double> coords1, std::vector<double> coords2 );

    double electronRepulsionPrimitive( Primitive * a, QuantumNumbers * ta, std::vector<double> A,
                                       Primitive * b, QuantumNumbers * tb, std::vector<double> B,
                                       Primitive * c, QuantumNumbers * tc, std::vector<double> C,
                                       Primitive * d, QuantumNumbers * td, std::vector<double> D );

    double electronRepulsionCGO( ContractedGaussianOrbital * a, std::vector<double> A,
                                 ContractedGaussianOrbital * b, std::vector<double> B,
                                 ContractedGaussianOrbital * c, std::vector<double> C,
                                 ContractedGaussianOrbital * d, std::vector<double> D );

    void fillOverlapMatrix( );
    void fillKineticEnergyMatrix( );
    void fillNuclearAttractionMatrix( );
    void fillElectronRepulsionMatrix( );

    void saveOverlapMatrix( std::string filename = "" );
    void saveKineticEnergyMatrix( std::string filename = "" );
    void saveNuclearAttractionMatrix( std::string filename = "" );

    void showElectronAttractionMatrix();

    Atom * getAtom( int n ) { return &atoms.at(n); }

    void setBasis( Basis * basis ) { this->basis = basis; }
    void showAtoms();

private:
    std::string scratch = "scratch/";

    const int MAXLINE = 100;

    std::vector<Atom> atoms;
    Basis * basis;

    Eigen::MatrixXd overlapMatrix;
    Eigen::MatrixXd kineticEnergyMatrix;
    Eigen::MatrixXd nuclearAttractionMatrix;
    Eigen::Tensor<double, 4> electronRepulsionTensor;
};
