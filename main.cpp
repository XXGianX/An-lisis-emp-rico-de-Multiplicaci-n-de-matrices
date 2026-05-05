#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <string>
#include <algorithm>
#include "matriz.h"

using namespace std;
namespace fs = std::filesystem;

int main() {
    auto start_total = chrono::high_resolution_clock::now();
    vector<int> dimensiones = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    const int repeticiones = 15;
    const int umbral_config = 16;
    int n0_detectado = -1;

    fs::create_directory("data");
    fs::create_directory("figures");

    ofstream csvFile("data/tiempos.csv");
    csvFile << "n,algoritmo,tiempo_ms\n";

    for (int n : dimensiones) {
        vector<double> tiemposE, tiemposS;
        for (int r = 0; r <= repeticiones; ++r) {
            Matrix A = generarMatriz(n, 100 + r);
            Matrix B = generarMatriz(n, 200 + r);

            auto s1 = chrono::high_resolution_clock::now();
            multiplicarEstandar(A, B, n);
            auto e1 = chrono::high_resolution_clock::now();
            double tE = chrono::duration<double, milli>(e1 - s1).count();

            int n_p2 = proximaPotenciaDe2(n);
            Matrix Ap = aplicarPadding(A, n, n_p2);
            Matrix Bp = aplicarPadding(B, n, n_p2);
            auto s2 = chrono::high_resolution_clock::now();
            multiplicarStrassen(Ap, Bp, n_p2, umbral_config);
            auto e2 = chrono::high_resolution_clock::now();
            double tS = chrono::duration<double, milli>(e2 - s2).count();

            if (r > 0) {
                tiemposE.push_back(tE);
                tiemposS.push_back(tS);
                csvFile << n << ",Estandar," << tE << "\n";
                csvFile << n << ",Strassen," << tS << "\n";
            }
        }

        sort(tiemposE.begin(), tiemposE.end());
        sort(tiemposS.begin(), tiemposS.end());
        if (n0_detectado == -1 && tiemposS[repeticiones/2] < tiemposE[repeticiones/2]) {
            n0_detectado = n;
        }
    }
    csvFile.close();

    auto end_total = chrono::high_resolution_clock::now();
    double total_seg = chrono::duration<double>(end_total - start_total).count();

    char marcaN0[256];
    if (n0_detectado != -1) {
        snprintf(marcaN0, sizeof(marcaN0),
                 "set arrow from %d,graph 0 to %d,graph 1 nohead lc rgb 'red' lw 2 dt 2\n"
                 "set label 'N0=%d' at %d,graph 0.1 tc rgb 'red' offset 1,0\n",
                 n0_detectado, n0_detectado, n0_detectado, n0_detectado);
    } else {
        snprintf(marcaN0, sizeof(marcaN0),
                 "set label 'N0 no detectado' at graph 0.05,0.93 tc rgb 'red'\n");
    }

    char script[2048];
    snprintf(script, sizeof(script),
        "set terminal pngcairo size 1000,620 enhanced font 'Monospace,11'\n"
        "set output 'figures/analisis_n0.png'\n"
        "set datafile separator ','\n"
        "set grid\n"
        "set logscale xy 2\n"
        "set title 'Estandar O(n^3) vs Strassen O(n^{2.81}) | umbral=%d | %.2fs'\n"
        "set xlabel 'n'\n"
        "set ylabel 'Tiempo (ms)'\n"
        "%s"
        "plot '< grep Estandar data/tiempos.csv' u 1:3 smooth unique "
             "with linespoints lc rgb '#185FA5' lw 2.5 pt 5 title 'Estandar',\\\n"
        "     '< grep Strassen data/tiempos.csv' u 1:3 smooth unique "
             "with linespoints lc rgb '#D85A30' lw 2.5 pt 7 title 'Strassen'\n",
        umbral_config, total_seg, marcaN0);

    //Ejecutar Gnuplot mediante popen para guardar y luego abrir ventana
    FILE* gp = popen("gnuplot", "w");
    if (!gp) { 
        fprintf(stderr, "Error: gnuplot no encontrado\n"); 
        return 1; 
    }
    fputs(script, gp);
    pclose(gp);

    // Abrir la ventana 
    system("gnuplot -p -e \"set datafile separator ','; set grid; set logscale xy 2; "
           "set title 'Analisis Finalizado'; "
           "plot '< grep Estandar data/tiempos.csv' u 1:3 smooth unique w lp title 'Estandar', "
           "'< grep Strassen data/tiempos.csv' u 1:3 smooth unique w lp title 'Strassen'\"");

    return 0;
}