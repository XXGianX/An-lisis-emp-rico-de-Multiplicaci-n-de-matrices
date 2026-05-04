#include "matriz.h"
#include <random>
#include <cmath>

using namespace std;

// Genera una matriz con valores aleatorios basados en una semilla
Matrix generarMatriz(int n, unsigned int seed) {
    Matrix m(n, vector<double>(n));
    mt19937 gen(seed); 
    uniform_real_distribution<> dis(1.0, 5.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            m[i][j] = dis(gen);
    return m;
}

// Extrae una submatriz de tamaño n x n desde una posición específica
void copiarSubmatriz(const Matrix& origen, Matrix& destino, int filaInicio, int colInicio, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            destino[i][j] = origen[filaInicio + i][colInicio + j];
        }
    }
}

// Inserta una submatriz en una matriz de mayor tamaño
void pegarSubmatriz(const Matrix& origen, Matrix& destino, int filaInicio, int colInicio, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            destino[filaInicio + i][colInicio + j] = origen[i][j];
        }
    }
}

// Calcula la dimensión mínima 2^k necesaria para el algoritmo de Strassen
int proximaPotenciaDe2(int n) {
    return pow(2, ceil(log2(n)));
}

// Incrementa el tamaño de la matriz rellenando con ceros
Matrix aplicarPadding(const Matrix& M, int n, int n_nuevo) {
    Matrix res(n_nuevo, vector<double>(n_nuevo, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            res[i][j] = M[i][j];
    return res;
}

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

// Multiplicación estándar con orden i-k-j para optimizar el uso de caché
Matrix multiplicarEstandar(const Matrix& A, const Matrix& B, int n) {
    Matrix C(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Implementación recursiva del algoritmo de Strassen
Matrix multiplicarStrassen(const Matrix& A, const Matrix& B, int n, int umbral) {
    if (n <= umbral) return multiplicarEstandar(A, B, n);

    int k = n / 2;
    Matrix A11(k, vector<double>(k)), A12(k, vector<double>(k)), A21(k, vector<double>(k)), A22(k, vector<double>(k));
    Matrix B11(k, vector<double>(k)), B12(k, vector<double>(k)), B21(k, vector<double>(k)), B22(k, vector<double>(k));

    copiarSubmatriz(A, A11, 0, 0, k); copiarSubmatriz(A, A12, 0, k, k);
    copiarSubmatriz(A, A21, k, 0, k); copiarSubmatriz(A, A22, k, k, k);
    copiarSubmatriz(B, B11, 0, 0, k); copiarSubmatriz(B, B12, 0, k, k);
    copiarSubmatriz(B, B21, k, 0, k); copiarSubmatriz(B, B22, k, k, k);

    // Productos intermedios P1 - P7
    Matrix P1 = multiplicarStrassen(sumar(A11, A22, k), sumar(B11, B22, k), k, umbral);
    Matrix P2 = multiplicarStrassen(sumar(A21, A22, k), B11, k, umbral);
    Matrix P3 = multiplicarStrassen(A11, restar(B12, B22, k), k, umbral);
    Matrix P4 = multiplicarStrassen(A22, restar(B21, B11, k), k, umbral);
    Matrix P5 = multiplicarStrassen(sumar(A11, A12, k), B22, k, umbral);
    Matrix P6 = multiplicarStrassen(restar(A21, A11, k), sumar(B11, B12, k), k, umbral);
    Matrix P7 = multiplicarStrassen(restar(A12, A22, k), sumar(B21, B22, k), k, umbral);

    // Construcción de cuadrantes de la matriz resultante
    Matrix C(n, vector<double>(n));
    Matrix C11 = sumar(restar(sumar(P1, P4, k), P5, k), P7, k);
    Matrix C12 = sumar(P3, P5, k);
    Matrix C21 = sumar(P2, P4, k);
    Matrix C22 = sumar(restar(sumar(P1, P3, k), P2, k), P6, k);

    pegarSubmatriz(C11, C, 0, 0, k); pegarSubmatriz(C12, C, 0, k, k);
    pegarSubmatriz(C21, C, k, 0, k); pegarSubmatriz(C22, C, k, k, k);

    return C;
}
