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

	int i;
	int j;
	int k;
};
