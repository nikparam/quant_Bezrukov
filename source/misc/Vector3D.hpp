#pragma once

#include <iostream>
#include <vector>

class Vector3D
{
public:
	Vector3D( double x, double y, double z );
	
	Vector3D & Vector3D::operator=( const Vector3D & other);
	Vector3D & Vector3D::operator+=( const Vector3D & other );

private:
	std::vector<double> coords;
};

