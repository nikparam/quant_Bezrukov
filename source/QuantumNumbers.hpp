#pragma once

#include <vector>

using std::vector;

class QuantumNumbers
{
public:
	explicit QuantumNumbers()
	{
	}

	QuantumNumbers( int i, int j, int k ) :
		i(i), j(j), k(k)
	{
	}

	int getAngularMomentum() const
	{
		return i + j + k;
	}

    QuantumNumbers change_i( const int v )
    {
        return QuantumNumbers(i + v, j, k);
    }

    QuantumNumbers change_j( const int v )
    {
        return QuantumNumbers(i, j + v, k);
    }

    QuantumNumbers change_k( const int v )
    {
        return QuantumNumbers(i, j, k + v);
    }

	int i;
	int j;
	int k;
};
