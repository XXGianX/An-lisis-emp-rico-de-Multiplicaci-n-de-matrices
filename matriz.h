#ifndef MATRIZ_H
#define MATRIZ_H

#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix generarMatriz(int n, unsigned int seed);
bool   sonIguales(const Matrix& A, const Matrix& B, int n);

void copiarSubmatriz(const Matrix& src, Matrix& dst, int fila, int col, int n);
void pegarSubmatriz (const Matrix& src, Matrix& dst, int fila, int col, int n);

int    proximaPotenciaDe2(int n);
Matrix aplicarPadding(const Matrix& M, int n, int nNuevo);
Matrix sumar (const Matrix& A, const Matrix& B, int n);
Matrix restar(const Matrix& A, const Matrix& B, int n);

Matrix multiplicarEstandar(const Matrix& A, const Matrix& B, int n);
Matrix multiplicarStrassen(const Matrix& A, const Matrix& B, int n, int umbral);

#endif