#include <iostream>
using namespace std;

#ifndef __Main_cpp__
#define __Main_cpp__

int main() {
	double* a = new double[100];
	memset(a, 1, sizeof(double) * 100);
	for (int i = 0; i < 100; i++) {
		std::cout << a[i] << " ";
	}
	delete[] a;
	return 0;
}

#endif