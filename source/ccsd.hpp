#pragma once

#include <iostream>
#include <cmath>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>

#include "ccsd_utilities.hpp"
#include "./molecule.hpp"

class CCSD
{
public:
    CCSD( int size_, int nocc, int nvirt, CCSD_Utilities& utilities) : size_(size_), nocc(nocc), nvirt(nvirt), utilities(utilities)
    {
        std::cout << "(CCSD) size_: " << size_ << "; nocc: " << nocc << "; nvirt: " << nvirt << std::endl;
    }

    void initialize();

    void fill_Dia();
    void fill_Dijab();

    void fill_initial_t1();
    void fill_initial_t2( Eigen::VectorXd const & HF_OrbitalEnergies );
    double test_MP2_Energy();

    void fill_tau();
    void fill_tilda_tau();

    void fill_Fae();
    void fill_Fme();
    void fill_Fmi();

    void fill_Wmnij();
    void fill_Wabef();
    void fill_Wmbej();

    void update_t1();
    void update_t2();

    double computeCCSD_correction();

    void preparation( Molecule const & molecule );
    void iterate();
    double run();
    double run_diis();

    double find_max( Eigen::MatrixXd const & B );

    Eigen::MatrixXd const & get_Dia() { return Dia; }
    Eigen::Tensor<double, 4> const & get_Dijab() { return Dijab; }

    Eigen::Tensor<double, 4> const & get_tau() { return tau; }
    Eigen::Tensor<double, 4> const & get_tilda_tau() { return tilda_tau; }

    Eigen::MatrixXd const & get_Fae() { return Fae; }
    Eigen::MatrixXd const & get_Fme() { return Fme; }
    Eigen::MatrixXd const & get_Fmi() { return Fmi; }

    Eigen::Tensor<double, 4> const & get_Wmnij() { return Wmnij; }
    Eigen::Tensor<double, 4> const & get_Wabef() { return Wabef; }
    Eigen::Tensor<double, 4> const & get_Wmbej() { return Wmbej; }

    Eigen::MatrixXd const & get_t1_updated() { return t1_updated; }
    Eigen::Tensor<double, 4> const & get_t2_updated() { return t2_updated; }

    inline int delta(int i, int j) { return (i == j); }

private:
    int size_;
    int nocc;
    int nvirt;

    CCSD_Utilities& utilities;

    // максимальное количество итераций CCSD
    const int maxiter = 40;
    // сходимость по энергии для CCSDcorr_E
    const double E_conv = 1.0e-12;
    // сходимость по энергии для CCSDcorr_E с DIIS
    const double E_conv_diis = 1.0e-12;
    // максимальное количество последних итераций, учитываемых в DIIS
    const int max_diis = 8;

    Eigen::MatrixXd Dia;
    Eigen::Tensor<double, 4> Dijab;

    Eigen::MatrixXd t1;
    Eigen::Tensor<double, 4> t2;

    Eigen::Tensor<double, 4> tau;
    Eigen::Tensor<double, 4> tilda_tau;

    Eigen::Tensor<double, 4> Wmnij;
    Eigen::Tensor<double, 4> Wabef;
    Eigen::Tensor<double, 4> Wmbej;

    Eigen::MatrixXd Fae;
    Eigen::MatrixXd Fme;
    Eigen::MatrixXd Fmi;

    Eigen::MatrixXd t1_updated;
    Eigen::Tensor<double, 4> t2_updated;
};

