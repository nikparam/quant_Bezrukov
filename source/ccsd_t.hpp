#pragma once

#include <iostream>
#include <cmath>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>

#include "ccsd_utilities.hpp"

class CCSD_T
{
public:
    CCSD_T( int size_, int nocc, int nvirt, CCSD_Utilities& utilities ) : size_(size_), nocc(nocc), nvirt(nvirt), utilities(utilities)
    {
        std::cout << "(CCSD_T) size_: " << size_ << "; nocc: " << nocc << "; nvirt: " << nvirt << std::endl;
    }

    void initialize();

    void set_t1( Eigen::MatrixXd const & t1 ) { this->t1 = t1; }
    void set_t2( Eigen::Tensor<double, 4> const & t2 ) { this->t2 = t2; }

    void build_disconnected_triples();
    void build_connected_triples();
    void build_Dijkabc();

    double compute_perturbation();

    Eigen::Tensor<double, 6> & get_t3d() { return t3d; }
    Eigen::Tensor<double, 6> & get_t3c() { return t3c; }

private:
    int size_;
    int nocc;
    int nvirt;
    CCSD_Utilities& utilities;

    Eigen::MatrixXd t1;
    Eigen::Tensor<double, 4> t2;

    Eigen::Tensor<double, 6> t3d; // disconnected triples
    Eigen::Tensor<double, 6> t3c; // connected triples
    Eigen::Tensor<double, 6> Dijkabc;
};

