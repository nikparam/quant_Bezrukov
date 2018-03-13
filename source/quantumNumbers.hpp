#pragma once

#include <vector>

using std::vector;

class QuantumNumbers
{
public:
	// пустой конструктор
	explicit QuantumNumbers()
	{
	}

	QuantumNumbers( int i, int j, int k ) :
		i(i), j(j), k(k)
	{
	}

	int getAngularMomentum()
	{
		return i + j + k;
	}

	int i;
	int j;
	int k;
};
