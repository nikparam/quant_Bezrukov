#include <iostream>
#include <vector>

using namespace std;

struct Foo
{
	~Foo() {} 
	double i = 0;
};
 
int main()
{
	Foo * foo = new Foo;
	vector<Foo*> foo_vector;
	foo_vector.push_back( foo );
	
	foo = new Foo;
	foo_vector.push_back( foo );

	for ( Foo * obj : foo_vector )
		delete obj;

	return 0;
}
