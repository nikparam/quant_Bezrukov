#include <iostream>
#include <string>
#include <vector>
#include <locale>
#include <sstream>
#include <memory>

using namespace std;

struct Vec3
{
    int x, y, z;
    Vec3() : x(0), y(0), z(0) { }
    Vec3(int x, int y, int z) :x(x), y(y), z(z) { }
    friend std::ostream& operator<<(std::ostream& os, Vec3& v) {
        return os << '{' << "x:" << v.x << " y:" << v.y << " z:" << v.z  << '}';
    }
};
 
int main()
{
    // Use the default constructor.
    std::unique_ptr<Vec3> v1 = std::make_unique<Vec3>();

	vector<unique_ptr<Vec3>> vec;
	vec.push_back( std::move(v1) );
	

	return 0;
}
