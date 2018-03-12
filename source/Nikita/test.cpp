#include <iostream>

int main(){
	int a = 1, b = 2;
	int * p1 = &a;
	int * p2 = &b;

	std::cout << *p1 + 1 << " " << *p2 * 2 << std::endl;

	return 0;
}
