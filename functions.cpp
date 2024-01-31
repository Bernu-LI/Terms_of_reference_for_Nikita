#include "functions.h"

#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <set>
#include <complex>
#include <string>

void SoLAE(std::vector <std::vector <double>>& matrix) {
    std::mt19937 rd(std::time(nullptr));
    std::uniform_int_distribution <int> dist(-9, 9);
    for (auto& row : matrix) {
        for (auto& value : row) {
            value = dist(rd);
        }
    }

    for (auto& row : matrix) {
        for (int i = 0; i < 3; i++) {
            if (row[i] < 0) {
                std::cout << "(" << (int)row[i] << ")x" << i + 1 << " + ";
            }
            else if (row[i] > 0) {
                std::cout << (int)row[i] << "x" << i + 1 << " + ";
            }
        }
        if (row[3] < 0) {
            std::cout << "(" << (int)row[3] << ")x4";
        }
        else if (row[3] > 0) {
            std::cout << (int)row[3] << "x4";
        }
        std::cout << " = " << (int)row[4] << "\n";
    }
}

std::vector <double> SoLAE_solution(std::vector <std::vector <double>>& matrix) {
    int numRows = 3;
    int numColumns = 5;

    // Forward elimination
    for (int p = 0; p < numRows; p++) {
        // Make sure pivot is non-zero
        if (std::abs(matrix[p][p]) < 1e-10) {
            std::cerr << "Mathematical error!" << std::endl;
            exit(1); // For this example, we'll just exit if we have a singular matrix.
        }

        // Eliminate the p-th column
        for (int i = p + 1; i < numRows; i++) {
            int alpha = matrix[i][p] / matrix[p][p];
            for (int j = p; j < numColumns; j++) {
                matrix[i][j] -= alpha * matrix[p][j];
            }
        }
    }

    // Backward substitution
    std::vector<double> answer(numRows);
    for (int i = numRows - 1; i >= 0; i--) {
        double sum = 0;
        for (size_t j = i + 1; j < numRows; j++) {
            sum += matrix[i][j] * answer[j];
        }
        answer[i] = (matrix[i][numColumns - 1] - sum) / matrix[i][i];
    }

    return answer;
}

void determinant4x4(std::vector <std::vector <int>>& det) {
    std::set <int> posOfZeroEls;
    std::mt19937 rd(std::time(nullptr));
    std::uniform_int_distribution dist(0, 15);
    while (posOfZeroEls.size() < 8) {
        int position = dist(rd);
        posOfZeroEls.insert(position);
    }

    std::uniform_int_distribution dist1(-9, 9);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int position = i * 4 + j;
            if (posOfZeroEls.find(position) == posOfZeroEls.end()) {
                det[i][j] = dist1(rd);
            }
        }
    }

    for (auto& row : det) {
        std::cout << "| ";
        for (auto& el : row) {
            if (el >= 0) {
                std::cout << " " << el << " |";
            }
            else {
                std::cout << el << " |";
            }
        }
        std::cout << "\n";
    }
}

int det3x3_calculation(std::vector <std::vector <int>>& det) {
    return (det[0][0] * det[1][1] * det[2][2] + det[0][2] * det[1][0] * det[2][1] + det[2][0] * det[0][1] * det[1][2]) -
            (det[0][2] * det[1][1] * det[2][0] + det[0][0] * det[2][1] * det[1][2] + det[2][2] * det[1][0] * det[0][1]);
}

int det4x4_calculation(std::vector <std::vector <int>>& det) {
    int calculation = 0;
    for (int p = 0; p < 4; p++) {
        std::vector <std::vector <int>> tempDet(3, std::vector <int>(3));
        for (int i = 1; i < 4; i++) {
            int columnCounter = 0;
            for (int j = 0; j < 4; j++) {
                if (j != p) {
                    tempDet[i - 1][columnCounter] = det[i][j];
                    columnCounter++;
                }
            }
        }
        int minorDeterminant = det3x3_calculation(tempDet);
        if (p % 2 == 0) {
            calculation += det[0][p] * minorDeterminant;
        }
        else {
            calculation -= det[0][p] * minorDeterminant;
        }
    }
    return calculation;
}

void determinantnxn(std::vector <std::vector <int>>& det) {
    std::mt19937 rd(std::time(nullptr));
    std::uniform_int_distribution vDist(-9, 9);
    for (auto& row : det) {
        for (auto& value : row) {
            value = vDist(rd);
        }
    }

    for (auto& row : det) {
        std::cout << "| ";
        for (auto& value : row) {
            if (value < 0) {
                std::cout << value << "| ";
            }
            else {
                std::cout << value << " | ";
            }
        }
        std::cout << "\n";
    }
}

int detnxn_calculation(std::vector <std::vector <int>>& det) {
    int n = det.size(); // Size of the matrix (n x n)
    int calculation = 0;

    // Base case: when matrix is 1x1
    if (n == 1) {
        return det[0][0];
    }

    // Base case: when matrix is 2x2
    if (n == 2) {
        return det[0][0] * det[1][1] - det[0][1] * det[1][0];
    }

    std::vector <std::vector <int>> subdet(n - 1, std::vector <int>(n - 1));
    // Recursive case: expand along the first row
    for (int x = 0; x < n; x++) {
        int subi = 0; // submatrix's i value

        // Construct the submatrix
        for (int i = 1; i < n; i++) {
            int subj = 0; // submatrix's j value

            for (int j = 0; j < n; j++) {
                if (j == x) {
                    continue; // skip column of the current element
                }

                subdet[subi][subj] = det[i][j];
                subj++;
            }

            subi++;
        }

        // Recursive call
        int subDet = detnxn_calculation(subdet);

        // Alternating signs for the expansion
        if (x % 2 == 0) {
            calculation += det[0][x] * subDet;
        } else {
            calculation -= det[0][x] * subDet;
        }
    }

    return calculation;
}

void complexNumbers() {
    std::vector <int> coefficients(6);
    std::mt19937 rd(std::time(nullptr));
    std::uniform_int_distribution dist(-9, 9);
    for (auto& value : coefficients) {
        while (value == 0) {
            value = dist(rd);
        }
    }
    std::cout << "conj("; coefficients[0] > 0 ? std::cout << coefficients[0] : std::cout << "(" << coefficients[0] << ")"; std::cout << " + "; coefficients[1] > 0 ? std::cout << coefficients[1] : std::cout << "(" << coefficients[1] << ")"; std::cout << "i)";
    std::cout << " * ";
    std::cout << "("; coefficients[2] > 0 ? std::cout << coefficients[2] : std::cout << "(" << coefficients[2] << ")"; std::cout << " + "; coefficients[3] > 0 ? std::cout << coefficients[3] : std::cout << "(" << coefficients[3] << ")"; std::cout << "i)";
    std::cout << " / ";
    std::cout << "("; coefficients[4] > 0 ? std::cout << coefficients[4] : std::cout << "(" << coefficients[4] << ")"; std::cout << " + "; coefficients[5] > 0 ? std::cout << coefficients[5] : std::cout << "(" << coefficients[5] << ")"; std::cout << "i)";
    std::cout << "\n";

    std::complex <double> z1(coefficients[0], coefficients[1]);
    std::complex <double> z2(coefficients[2], coefficients[3]);
    std::complex <double> z3(coefficients[4], coefficients[5]);
    std::complex <double> result = std::conj(z1) * z2 / z3;
    std::cout << result.real() << " + " << (result.imag() < 0 ? "(" + std::to_string(result.imag()) + ")" : std::to_string(result.imag())) << "i";
}

void call() {
    std::cout << "1.\n";
    std::vector <std::vector <double>> matrix(3, std::vector <double>(5));
    SoLAE(matrix);
    std::vector <double> answer = SoLAE_solution(matrix);
    for (int i = 0; i < 4; i++) {
        std::cout << "x" << i + 1 << " = " << answer[i] << "\n";
    }
    std::cout << "\n";

    std::cout << "3.\n";
    std::vector <std::vector <int>> det(4, std::vector <int>(4));
    determinant4x4(det);
    std::cout << det4x4_calculation(det) ;
    std::cout << "\n\n";

    std::cout << "4.\n";
    std::mt19937 rd(std::time(nullptr));
    std::uniform_int_distribution dist(1, 10);
    int n = dist(rd);
    std::vector <std::vector <int>> det1(n, std::vector <int>(n));
    determinantnxn(det1);
    std::cout << detnxn_calculation(det1);
    std::cout << "\n\n";

    std::cout << "7.\n";
    complexNumbers();
}
