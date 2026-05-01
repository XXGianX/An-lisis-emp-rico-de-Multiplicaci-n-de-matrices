#ifndef MATRIZ_H
#define MATRIZ_H

#include <vector>

typedef std::vector<std::vector<double>> Matrix;

Matrix generarMatriz(int n);
Matrix multiplicarEstandar(const Matrix& A, const Matrix& B, int n);

#endif
