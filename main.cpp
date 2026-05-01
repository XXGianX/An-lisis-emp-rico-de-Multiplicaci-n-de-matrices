#include <iostream>
#include <chrono>
#include <fstream>
#include <numeric>
#include "matriz.h"

using namespace std; 

int main() {
    ofstream dataFile("data/tiempos_estandar.dat");
    int repeticiones = 5;

    // Evaluamos tamaños n (potencias de 2)
    for (int n = 2; n <= 512; n *= 2) {
        vector<double> tiempos;
        
        for (int r = 0; r < repeticiones; ++r) {
            Matrix A = generarMatriz(n);
            Matrix B = generarMatriz(n);
            
            auto start = chrono::high_resolution_clock::now();
            multiplicarEstandar(A, B, n);
            auto end = chrono::high_resolution_clock::now();
            
            chrono::duration<double> diff = end - start;
            tiempos.push_back(diff.count());
        }
        
        double promedio = accumulate(tiempos.begin(), tiempos.end(), 0.0) / repeticiones;
        
        dataFile << n << " " << promedio << "\n";
        cout << n << "\t" << promedio << "s" << endl;
    }
    dataFile.close();

    // Comando gnuplot optimizado para limpieza visual
    system("gnuplot -p -e \"set title 'Rendimiento Algoritmo Estándar'; \
                         set xlabel 'Tamaño de Matriz (n)'; \
                         set ylabel 'Tiempo Promedio (s)'; \
                         set grid; \
                         set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5; \
                         plot 'data/tiempos_estandar.dat' with linespoints ls 1 title 'O(n^3)'\"");

    return 0;
}
