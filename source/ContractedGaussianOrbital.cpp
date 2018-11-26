#include "ContractedGaussianOrbital.hpp"

void ContractedGaussianOrbital::add_primitive( double alpha, double exponent )
{
	primitives.emplace_back( alpha, exponent );
	primitives.end()[-1].normalize( qNumbers );
}

void ContractedGaussianOrbital::normalize()
{
    int L = qNumbers.getAngularMomentum();
    double prefactor = std::pow(M_PI, 1.5) * MathUtils::doubleFactorial(qNumbers.i - 1) * MathUtils::doubleFactorial(qNumbers.j - 1) * MathUtils::doubleFactorial(qNumbers.k - 1) / std::pow(2, L);

    double N = 0.0;

    size_t size_ = primitives.size();
    for ( size_t ia = 0; ia < size_; ++ia )
    {
        for ( size_t ib = 0; ib < size_; ++ib )
        {
            double norm_a = primitives[ia].get_norm( qNumbers );
            double norm_b = primitives[ib].get_norm( qNumbers );
            double exp_a = primitives[ia].get_exponent();
            double exp_b = primitives[ib].get_exponent();

            N += norm_a * norm_b * exp_a * exp_b / std::pow(exp_a + exp_b, L + 1.5);
        }
    }

    for ( size_t i = 0; i < size_; ++i )
        std::cout << "i: " << i << "; norm: " << primitives[i].get_norm(qNumbers) << std::endl;

    N *= prefactor;
    N = std::pow(N, -0.5);

    std::cout << "(CGO normalize) N: " << N << std::endl;

    for ( size_t i = 0; i < size_; ++i )
        primitives[i].multiply_exponent( N );
}

void ContractedGaussianOrbital::show()
{
    std::cout << "(CGO)" << primitives.size() << std::endl;
    std::cout << "(CGO) QuantumNumbers: " << qNumbers.i << " " << qNumbers.j << " " << qNumbers.k << std::endl;

	for ( auto p : primitives )
		p.show();
}
