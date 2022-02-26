#include "ArrayIO.h"

bool BinaryDataIO::Write_Scalar_Array(const std::string& file_name, const Array<bool>& arr) {
	std::ofstream output(file_name, std::ios::binary);
	if (!output) return false;
	std::uint32_t n = (std::uint32_t)arr.size();
	File::Write_Binary<std::uint32_t>(output, n);
	bool* data = new bool[n];
#pragma omp parallel for
	for (int i = 0; i < (int)n; i++) data[i] = arr[i];
	File::Write_Binary_Array<bool>(output, data, n);
	delete[] data;
	output.close();
	return true;
}

bool BinaryDataIO::Read_Scalar_Array(const std::string& file_name, Array<bool>& arr)
{
	std::ifstream input(file_name, std::ios::binary);
	if (!input) return false;
	std::uint32_t n;
	File::Read_Binary<std::uint32_t>(input, n);
	arr.resize(n);
	bool* data = new bool[n];
	File::Read_Binary_Array<bool>(input, data, (int)n);
	input.close();
#pragma omp parallel for
	for (int i = 0; i < (int)n; i++) arr[i] = data[i];
	delete[] data;
	return true;
}