#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

class MP2
{
public:
    MP2( int size_, int occ, Eigen::Tensor<double, 4> const& electronRepulsionTensor, Eigen::MatrixXd const& matrixC,
         Eigen::VectorXd const& HF_OrbitalEnergies ) : size_(size_), occ(occ),
        electronRepulsionTensor(electronRepulsionTensor), matrixC(matrixC), HF_OrbitalEnergies(HF_OrbitalEnergies)
    {
    }

    void fillTwoElectronMOIntegrals();
    void fillTwoElectronMOIntegrals_eff();
    double computeMP2_correction( );

    Eigen::Tensor<double, 4> const& get_two_electron_MO_integrals() const { return twoElectronMOIntegrals; }

private:
    int size_;
    int occ;

    Eigen::Tensor<double, 4> twoElectronMOIntegrals;

    Eigen::Tensor<double, 4> const& electronRepulsionTensor;
    Eigen::MatrixXd const& matrixC;
    Eigen::VectorXd const& HF_OrbitalEnergies;
};
