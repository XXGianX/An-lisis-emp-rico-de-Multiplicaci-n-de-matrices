#ifndef MATRIZ_H
#define MATRIZ_H

#include <vector>

typedef std::vector<std::vector<double>> Matrix;

// Utilidades
Matrix generarMatriz(int n, unsigned int seed);
bool sonIguales(const Matrix& A, const Matrix& B, int n);

// Algoritmos
Matrix multiplicarEstandar(const Matrix& A, const Matrix& B, int n);
Matrix multiplicarStrassen(const Matrix& A, const Matrix& B, int n);

#endif
