#include <iostream>
#include <chrono>
#include <fstream>
#include <numeric>
#include "matriz.h"

using namespace std;

int main() {
    ofstream dataFile("data/comparativa.dat");
    int repeticiones = 5;
    unsigned int semillaBase = 42;

    cout << "--- Benchmark: Estándar vs Strassen (Ryzen 5 5600H) ---" << endl;
    cout << "N\tEstándar(s)\tStrassen(s)\tEstado" << endl;
    cout << "--------------------------------------------" << endl;

    for (int n = 2; n <= 128; n *= 2) {
        vector<double> tStd, tStr;
        bool correcto = true;

        for (int r = 0; r < repeticiones; ++r) {
            Matrix A = generarMatriz(n, semillaBase + r);
            Matrix B = generarMatriz(n, semillaBase + r + 100);

            auto s1 = chrono::high_resolution_clock::now();
            Matrix resStd = multiplicarEstandar(A, B, n);
            auto e1 = chrono::high_resolution_clock::now();
            tStd.push_back(chrono::duration<double>(e1 - s1).count());

            auto s2 = chrono::high_resolution_clock::now();
            Matrix resStr = multiplicarStrassen(A, B, n);
            auto e2 = chrono::high_resolution_clock::now();
            tStr.push_back(chrono::duration<double>(e2 - s2).count());

            if (!sonIguales(resStd, resStr, n)) correcto = false;
        }

        double pStd = accumulate(tStd.begin(), tStd.end(), 0.0) / repeticiones;
        double pStr = accumulate(tStr.begin(), tStr.end(), 0.0) / repeticiones;

        cout << n << "\t" << pStd << "\t" << pStr << "\t" << (correcto ? "OK" : "FAIL") << endl;
        dataFile << n << " " << pStd << " " << pStr << "\n";
    }
    dataFile.close();

    cout << "\nGenerando gráfica..." << endl;
    // Capturamos el retorno de system para silenciar al compilador
    int res = system("gnuplot -p -e \"set title 'Estándar vs Strassen'; \
                      set xlabel 'n'; set ylabel 'Tiempo (s)'; set grid; \
                      plot 'data/comparativa.dat' using 1:2 with linespoints title 'Estándar', \
                           'data/comparativa.dat' using 1:3 with linespoints title 'Strassen'\"");
    (void)res; 

    return 0;
}
