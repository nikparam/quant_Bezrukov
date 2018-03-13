#include <iostream>
#include <vector>

#include "MathUtils.hpp"

using namespace std;

struct Foo
{
	~Foo() {} 
	double i = 0;
};
 
int main()
{
	cout << "1!!: " << MathUtils::doubleFactorial( 1 ) << endl;
	cout << "2!!: " << MathUtils::doubleFactorial( 2 ) << endl;
	cout << "3!!: " << MathUtils::doubleFactorial( 3 ) << endl;
	cout << "4!!: " << MathUtils::doubleFactorial( 4 ) << endl;
	cout << "5!!: " << MathUtils::doubleFactorial( 5 ) << endl;
	cout << "6!!: " << MathUtils::doubleFactorial( 6 ) << endl;
	cout << "7!!: " << MathUtils::doubleFactorial( 7 ) << endl;
	
	return 0;
}
