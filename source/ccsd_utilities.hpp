#pragma once

// Будем наследовать CCSD и CCSD_T от этого класса
// Сюда поместим SO матрицу Фока, матрицу антисимметризованных двухэлектронных интегралов
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>

class CCSD_Utilities
{
public:
    CCSD_Utilities() {}

    CCSD_Utilities( int size_, int nocc, int nvirt ) : size_(size_), nocc(nocc), nvirt(nvirt)
    {
        std::cout << "(CCSD_Utilities) size_: " << size_ << "; nocc: " << nocc << "; nvirt: " << nvirt << std::endl;
        initialize_utilities();
    }

    void initialize_utilities();

    int size_;
    int nocc;
    int nvirt;

    void fillAS_MO_TwoElectronIntegrals( Eigen::Tensor<double, 4> const & twoElectronMOIntegrals );
    void fillSOHcore( Eigen::MatrixXd const & matrixC, Eigen::MatrixXd const & matrixHcore );
    void fillSOFock();

    Eigen::MatrixXd const & get_Hcore() { return SOHcoreMatrix; }
    Eigen::MatrixXd const & get_SOF() { return SOFockMatrix; }

    Eigen::Tensor<double, 4> ASTwoElectronMOIntegrals;
    Eigen::MatrixXd SOHcoreMatrix;
    Eigen::MatrixXd SOFockMatrix;
};
