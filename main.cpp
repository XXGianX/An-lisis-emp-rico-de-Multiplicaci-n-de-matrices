#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include "matriz.h"

using namespace std;

int main() {
    // Configuración de parámetros experimentales
    vector<int> dimensiones = {2, 4, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256};
    const int repeticiones = 15;
    const int umbral_config = 32;

    ofstream csvFile("data/tiempos.csv");
    if (!csvFile.is_open()) return 1;

    csvFile << "n,algoritmo,repeticion,tiempo_ms\n";

    for (int n : dimensiones) {
        cout << "Ejecutando pruebas para n = " << n << "..." << endl;
        
        for (int r = 1; r <= repeticiones; ++r) {
            Matrix A = generarMatriz(n, 100 + r);
            Matrix B = generarMatriz(n, 200 + r);

            // Medición de rendimiento: Algoritmo Estándar
            auto s1 = chrono::high_resolution_clock::now();
            multiplicarEstandar(A, B, n);
            auto e1 = chrono::high_resolution_clock::now();
            csvFile << n << ",Estandar," << r << "," 
                    << chrono::duration<double, milli>(e1 - s1).count() << "\n";

            // Medición de rendimiento: Algoritmo Strassen
            int n_p2 = proximaPotenciaDe2(n);
            Matrix Ap = aplicarPadding(A, n, n_p2);
            Matrix Bp = aplicarPadding(B, n, n_p2);

            auto s2 = chrono::high_resolution_clock::now();
            multiplicarStrassen(Ap, Bp, n_p2, umbral_config);
            auto e2 = chrono::high_resolution_clock::now();
            csvFile << n << ",Strassen," << r << "," 
                    << chrono::duration<double, milli>(e2 - s2).count() << "\n";
        }
    }

    csvFile.close();
    
    // Generación de reporte visual (Gnuplot)
    // Se captura el valor de retorno para cumplir con el atributo warn_unused_result
    // Reemplaza el bloque de gnuplot en tu main.cpp por este:
    int status = system("gnuplot -p -e \"set title 'Rendimiento: Estandar vs Strassen (Promediado)'; \
                    set datafile separator ','; set xlabel 'Dimension (n)'; \
                    set ylabel 'Tiempo (ms)'; set grid; \
                    plot 'data/tiempos.csv' u 1:4 smooth unique with linespoints title 'Estandar', \
                         'data/tiempos.csv' u 1:4 smooth unique with linespoints title 'Strassen'\"");
    return (status != -1) ? 0 : 1;
}
