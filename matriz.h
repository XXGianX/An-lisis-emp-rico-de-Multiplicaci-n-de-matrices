#ifndef MATRIZ_H
#define MATRIZ_H

#include <vector>

using namespace std;

typedef vector<vector<double>> Matrix;

Matrix generarMatriz(int n, unsigned int seed);
bool sonIguales(const Matrix& A, const Matrix& B, int n);

int proximaPotenciaDe2(int n);
Matrix aplicarPadding(const Matrix& M, int n, int n_nuevo);


Matrix sumar(const Matrix& A, const Matrix& B, int n);
Matrix restar(const Matrix& A, const Matrix& B, int n);


Matrix multiplicarEstandar(const Matrix& A, const Matrix& B, int n);
Matrix multiplicarStrassen(const Matrix& A, const Matrix& B, int n, int umbral);

#endif
