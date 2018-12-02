#pragma once

#include "atom.hpp"
#include "Basis.hpp"
#include "ContractedGaussianOrbital.hpp"
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/MatrixFunctions>

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

    void fillGMatrix( );
    void fillFMatrix( );

    double find_max( const Eigen::MatrixXd & B );
    double compute_E0();
    double compute_Enuc();

    // SCF procedures
    void SCF_initialize( );
    double SCF( );
    double SCF_DIIS( );

    void saveOverlapMatrix( std::string filename = "" );
    void saveKineticEnergyMatrix( std::string filename = "" );
    void saveNuclearAttractionMatrix( std::string filename = "" );
    void saveXMatrix( std::string filename = "" );
    void saveGMatrix( std::string filename = "" );
    void saveFMatrix( std::string filename = "" );

    void setCharge( const int icharge = 0 );
    void setOutput( std::string const & filename = "" );

    void showElectronRepulsionTensor();

    void makeInitialGuess( );

    Atom * getAtom( int n ) { return &atoms.at(n); }

    void setBasis( Basis * basis ) { this->basis = basis; }
    void showAtoms();

    int size( ) const;
    int delta(int i, int j) { return (i == j); }
    int get_charge() const { return charge; }

    Eigen::MatrixXd const & get_C() const { return matrixC; }
    Eigen::MatrixXd const & get_Hcore() const { return matrixHcore; }
    Eigen::VectorXd const & get_HF_OrbitalEnergies() const { return HF_OrbitalEnergies; }

    Eigen::Tensor<double, 4> const & get_two_electron_integrals() const { return electronRepulsionTensor; }

private:
    int charge;
    bool set_charge = false;

    std::string scratch = "scratch/";

    const int MAXLINE = 100;

    std::vector<Atom> atoms;
    Basis * basis;

    Eigen::MatrixXd overlapMatrix;
    Eigen::MatrixXd kineticEnergyMatrix;
    Eigen::MatrixXd nuclearAttractionMatrix;
    Eigen::Tensor<double, 4> electronRepulsionTensor;

    Eigen::MatrixXd matrixP;
    Eigen::MatrixXd matrixP_new; // получаемая на текущей итерации
    Eigen::MatrixXd matrixG;
    Eigen::MatrixXd matrixF;
    Eigen::MatrixXd matrixFprime;
    Eigen::MatrixXd matrixC;
    Eigen::MatrixXd matrixHcore;

    Eigen::VectorXd HF_OrbitalEnergies;

    std::ofstream outFile;

    const double E_conv = 1.0e-8; // сходимость SCF по энергии
    const double max_diis_size = 8;
    const int max_diis_iter = 15;
};
