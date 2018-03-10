#include <iostream>

int main(){
	int * p1 = nullptr;
	int * p2 = nullptr;
	int a, b;
	if ( std::cin >> a ) p1 = &a;
	if ( std::cin >> a ) p2 = &b;
	if ( * p1 != * p2 ) std::cout << "Success" << std::endl;
	return 0;
}
