#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "ccsd_utilities.hpp"

class CIS
{
public:
    CIS( int size_, int nocc, int nvirt, CCSD_Utilities & utilities) :
        size_(size_), nocc(nocc), nvirt(nvirt), utilities(utilities)
    {
        std::cout << "Initializing CIS class." << std::endl << "size_: " << size_ << "; nocc: " << nocc << "; nvirt: " << nvirt << std::endl;
    }

    void initialize();
    void fill_cis_matrix();
    Eigen::VectorXd diagonalize_cis_matrix();

    Eigen::MatrixXd & get_cis_matrix() { return cis_matrix; }

    inline int delta(int i, int j) { return i == j; }
private:
    int size_;
    int nocc;
    int nvirt;

    CCSD_Utilities & utilities;

    Eigen::MatrixXd cis_matrix;
};

