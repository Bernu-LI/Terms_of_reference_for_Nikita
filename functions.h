#ifndef functions
#define functions

#include <vector>

void SoLAE(std::vector <std::vector <double>>& matrix);
std::vector<double> SoLAE_solution (std::vector <std::vector <double>>& matrix);

void determinant4x4(std::vector <std::vector <int>>& det);
int det3x3_calculation(std::vector <std::vector <int>>& det);

void determinantnxn(std::vector <std::vector <int>>& det);
int detnxn_calculation(std::vector <std::vector <int>>& det);

void complexNumbers();

void call();

#endif