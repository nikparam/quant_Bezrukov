#pragma once

#include <iostream>
#include <vector>

using std::vector;

class Vector3D
{
public:
	Vector3D( double x, double y, double z )
	{
		coords.push_back( x );
		coords.push_back( y );
		coords.push_back( z );
	}

private:
	vector<double> coords;
};

//Vector3D & Vector3D::operator=( const Vector3D & other)
//{
	//if ( this != &other )
	//{
		//coords = other.coords;
		//return *this;
	//}

	//return *this;
//}

//Vector3D & Vector3D::operator+=( const Vector3D & other )
//{
	//for ( size_t i = 0; i < 3; i++ )
		//coords[i] += other.coords[i];
	
	//return *this;
//}

