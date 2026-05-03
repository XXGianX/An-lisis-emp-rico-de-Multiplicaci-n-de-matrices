#include "matriz.h"
#include <random>
#include <cmath>
#include <algorithm>

using namespace std;

// Generador con semilla para reproducibilidad
Matrix generarMatriz(int n, unsigned int seed) {
    Matrix m(n, vector<double>(n));
    mt19937 gen(seed); 
    uniform_real_distribution<> dis(1.0, 5.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m[i][j] = dis(gen);
    return m;
}

bool sonIguales(const Matrix& A, const Matrix& B, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (abs(A[i][j] - B[i][j]) > 1e-7) return false;
    return true;
}

// Auxiliares de Padding
int proximaPotenciaDe2(int n) {
    return pow(2, ceil(log2(n)));
}

Matrix aplicarPadding(const Matrix& M, int n, int n_nuevo) {
    Matrix res(n_nuevo, vector<double>(n_nuevo, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            res[i][j] = M[i][j];
    return res;
}

// Operaciones básicas
Matrix sumar(const Matrix& A, const Matrix& B, int n) {
    Matrix C(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

Matrix restar(const Matrix& A, const Matrix& B, int n) {
    Matrix C(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

// Algoritmo Estándar O(n^3)
Matrix multiplicarEstandar(const Matrix& A, const Matrix& B, int n) {
    Matrix C(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) { // Optimización de localidad: i-k-j
            for (int j = 0; j < n; ++j) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Algoritmo de Strassen con Umbral
Matrix multiplicarStrassen(const Matrix& A, const Matrix& B, int n, int umbral) {
    if (n <= umbral) { 
        return multiplicarEstandar(A, B, n);
    }

    int k = n / 2;
    Matrix A11(k, vector<double>(k)), A12(k, vector<double>(k)), A21(k, vector<double>(k)), A22(k, vector<double>(k));
    Matrix B11(k, vector<double>(k)), B12(k, vector<double>(k)), B21(k, vector<double>(k)), B22(k, vector<double>(k));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A11[i][j] = A[i][j];         A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];     A22[i][j] = A[i + k][j + k];
            B11[i][j] = B[i][j];         B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];     B22[i][j] = B[i + k][j + k];
        }
    }

    Matrix P1 = multiplicarStrassen(sumar(A11, A22, k), sumar(B11, B22, k), k, umbral);
    Matrix P2 = multiplicarStrassen(sumar(A21, A22, k), B11, k, umbral);
    Matrix P3 = multiplicarStrassen(A11, restar(B12, B22, k), k, umbral);
    Matrix P4 = multiplicarStrassen(A22, restar(B21, B11, k), k, umbral);
    Matrix P5 = multiplicarStrassen(sumar(A11, A12, k), B22, k, umbral);
    Matrix P6 = multiplicarStrassen(restar(A21, A11, k), sumar(B11, B12, k), k, umbral);
    Matrix P7 = multiplicarStrassen(restar(A12, A22, k), sumar(B21, B22, k), k, umbral);

    Matrix C(n, vector<double>(n));
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j] = P1[i][j] + P4[i][j] - P5[i][j] + P7[i][j];
            C[i][j + k] = P3[i][j] + P5[i][j];
            C[i + k][j] = P2[i][j] + P4[i][j];
            C[i + k][j + k] = P1[i][j] - P2[i][j] + P3[i][j] + P6[i][j];
        }
    }
    return C;
}
