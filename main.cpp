#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include "matriz.h"

using namespace std;

// Función para calcular la mediana y reducir ruido
double obtenerMediana(vector<double>& tiempos) {
    sort(tiempos.begin(), tiempos.end());
    int n = tiempos.size();
    if (n % 2 == 0) return (tiempos[n/2 - 1] + tiempos[n/2]) / 2.0;
    return tiempos[n/2];
}

int main() {
    ofstream dataFile("data/benchmark_final.dat");
    int repeticiones = 7; 
    int umbral_optimo = 64; // Basado en pruebas previas

    cout << "Benchmark para Gian (UNSA - Arequipa)" << endl;
    cout << "N\tEstándar(s)\tStrassen(s)\tResultado" << endl;
    cout << "------------------------------------------------" << endl;

    for (int n : {64, 128, 256, 512}) {
        vector<double> tStd, tStr;
        bool exito = true;

        for (int r = 0; r < repeticiones; ++r) {
            Matrix A_orig = generarMatriz(n, 42 + r);
            Matrix B_orig = generarMatriz(n, 99 + r);

            // Manejo de Padding
            int n_p2 = proximaPotenciaDe2(n);
            Matrix A_p = aplicarPadding(A_orig, n, n_p2);
            Matrix B_p = aplicarPadding(B_orig, n, n_p2);

            auto s1 = chrono::high_resolution_clock::now();
            Matrix resStd = multiplicarEstandar(A_orig, B_orig, n);
            auto e1 = chrono::high_resolution_clock::now();
            tStd.push_back(chrono::duration<double>(e1 - s1).count());

            auto s2 = chrono::high_resolution_clock::now();
            Matrix resStr_p = multiplicarStrassen(A_p, B_p, n_p2, umbral_optimo);
            auto e2 = chrono::high_resolution_clock::now();
            tStr.push_back(chrono::duration<double>(e2 - s2).count());

            // Verificar contra el estándar (solo la parte útil de la matriz)
            for(int i=0; i<n; i++)
                for(int j=0; j<n; j++)
                    if(abs(resStd[i][j] - resStr_p[i][j]) > 1e-7) exito = false;
        }

        double medStd = obtenerMediana(tStd);
        double medStr = obtenerMediana(tStr);

        cout << n << "\t" << medStd << "\t" << medStr << "\t" << (exito ? "CORRECTO" : "ERROR") << endl;
        dataFile << n << " " << medStd << " " << medStr << "\n";
    }

    dataFile.close();

    // Silenciando advertencias de Ubuntu/GCC
    cout << "\nGenerando gráfica comparativa..." << endl;
    int check = system("gnuplot -p -e \"set title 'Rendimiento: Estándar vs Strassen (Mediana)'; \
                        set xlabel 'n'; set ylabel 'Tiempo (s)'; set grid; \
                        plot 'data/benchmark_final.dat' using 1:2 with linespoints title 'Estándar', \
                             'data/benchmark_final.dat' using 1:3 with linespoints title 'Strassen'\"");
    (void)check; 

    return 0;
}
