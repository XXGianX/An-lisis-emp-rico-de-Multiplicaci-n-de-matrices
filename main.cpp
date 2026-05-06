#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <stdexcept>
#include "matriz.h"

using namespace std;
namespace fs = std::filesystem;

// Tipos
struct ResultadoN {
    int    n;
    double medianaEstandar, mediaEstandar;
    double medianaStrassen, mediaStrassen;
};

// Estadísticos
static double mediana(vector<double> v) {
    sort(v.begin(), v.end());
    int n = static_cast<int>(v.size());
    return (n % 2 == 0) ? (v[n/2 - 1] + v[n/2]) / 2.0 : v[n/2];
}

static double media(const vector<double>& v) {
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// Medición 

static ResultadoN ejecutarBenchmark(int n, int reps, int umbral, FILE* csv) {
    int  nPad    = proximaPotenciaDe2(n);
    bool padding = (nPad != n);

    // Warmup: descarta la primera ejecucion para estabilizar cache y branch predictor
    {
        Matrix A  = generarMatriz(n, 0), B = generarMatriz(n, 1);
        Matrix Ap = padding ? aplicarPadding(A, n, nPad) : A;
        Matrix Bp = padding ? aplicarPadding(B, n, nPad) : B;
        multiplicarEstandar(A,  B,  n);
        multiplicarStrassen(Ap, Bp, nPad, umbral);
    }

    vector<double> tE, tS;
    tE.reserve(reps);
    tS.reserve(reps);

    for (int r = 1; r <= reps; ++r) {
        Matrix A  = generarMatriz(n, 100 + r);
        Matrix B  = generarMatriz(n, 200 + r);
        Matrix Ap = padding ? aplicarPadding(A, n, nPad) : A;
        Matrix Bp = padding ? aplicarPadding(B, n, nPad) : B;

        auto s1 = chrono::high_resolution_clock::now();
        multiplicarEstandar(A, B, n);
        auto e1 = chrono::high_resolution_clock::now();

        auto s2 = chrono::high_resolution_clock::now();
        multiplicarStrassen(Ap, Bp, nPad, umbral);
        auto e2 = chrono::high_resolution_clock::now();

        double tEval = chrono::duration<double, milli>(e1 - s1).count();
        double tSval = chrono::duration<double, milli>(e2 - s2).count();

        tE.push_back(tEval);
        tS.push_back(tSval);

        fprintf(csv, "%d,Estandar,%.6f\n", n, tEval);
        fprintf(csv, "%d,Strassen,%.6f\n", n, tSval);
    }

    return { n, mediana(tE), media(tE), mediana(tS), media(tS) };
}

// Salida 

static void imprimirTabla(const vector<ResultadoN>& res, int umbral) {
    static const char SEP[] =
        "--------------------------------------------------------------------------\n";

    printf("\n  Benchmark — Estandar O(n^3) vs Strassen O(n^2.81)\n"
           "  Repeticiones: 15  |  Estadistico: mediana  |  Umbral: %d\n"
           "  AMD Ryzen 5 5600H · Ubuntu · g++ -O2\n\n", umbral);
    printf("%s", SEP);
    printf("%6s %16s %16s %16s %16s %9s\n",
           "n", "Est.med(ms)", "Est.avg(ms)", "Str.med(ms)", "Str.avg(ms)", "Ganador");
    printf("%s", SEP);

    for (const auto& r : res) {
        const char* ganador = (r.medianaStrassen < r.medianaEstandar) ? "Strassen" : "Estandar";
        printf("%6d %16.4f %16.4f %16.4f %16.4f %9s\n",
               r.n, r.medianaEstandar, r.mediaEstandar,
               r.medianaStrassen, r.mediaStrassen, ganador);
    }
    printf("%s", SEP);
}

static void escribirCSVMedianas(const vector<ResultadoN>& res, const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) throw runtime_error(string("No se pudo crear: ") + path);
    fprintf(f, "n,mediana_estandar,mediana_strassen\n");
    for (const auto& r : res)
        fprintf(f, "%d,%.6f,%.6f\n", r.n, r.medianaEstandar, r.medianaStrassen);
    fclose(f);
}

static int detectarN0(const vector<ResultadoN>& res) {
    // Busca el primer n donde Strassen supera consistentemente al estandar
    for (const auto& r : res)
        if (r.medianaStrassen < r.medianaEstandar) return r.n;
    return -1;
}

// Gnuplot

static void graficar(const char* csvMedianas, const char* outDir,
                     int n0, int umbral, double totalSeg) {
    char lineal[256], loglog[256], marcaN0[512], script[4096];

    snprintf(lineal, sizeof(lineal), "%s/comparacion_lineal.png", outDir);
    snprintf(loglog, sizeof(loglog), "%s/comparacion_loglog.png", outDir);

    if (n0 != -1) {
        snprintf(marcaN0, sizeof(marcaN0),
                 "set arrow from %d,graph 0 to %d,graph 1 nohead lc rgb '#D85A30' lw 2 dt 2\n"
                 "set label 'N0=%d' at %d,graph 0.08 tc rgb '#D85A30' "
                 "font 'Monospace Bold,10' offset 1,0\n",
                 n0, n0, n0, n0);
    } else {
        snprintf(marcaN0, sizeof(marcaN0),
                 "set label 'N0 no detectado en rango evaluado' "
                 "at graph 0.04,0.93 tc rgb '#D85A30' font 'Monospace,10'\n");
    }

    snprintf(script, sizeof(script),
        "set terminal pngcairo size 1000,620 enhanced font 'Monospace,11'\n"
        "set datafile separator ','\n"
        "set grid lc rgb '#dddddd' lt 1 lw 0.6\n"
        "set border lw 1.2\n"
        "set key top left box lw 0.8 samplen 2 spacing 1.3\n"
        "set pointsize 1.1\n"

        // Escala lineal
        "set output '%s'\n"
        "set title 'Estandar O(n^3) vs Strassen O(n^{2.81}) - escala lineal\\n"
            "umbral=%d - AMD Ryzen 5 5600H - Ubuntu - g++ -O2 - 15 reps - %.2fs' "
            "font 'Monospace Bold,11'\n"
        "set xlabel 'Dimension n' font 'Monospace,11'\n"
        "set ylabel 'Tiempo (ms) - mediana 15 muestras' font 'Monospace,11'\n"
        "set xtics rotate by -45\n"
        "%s"
        "plot '%s' using 1:2 with linespoints lc rgb '#185FA5' lw 2.5 pt 5 "
            "title 'Estandar O(n^3)',\\\n"
        "     '%s' using 1:3 with linespoints lc rgb '#D85A30' lw 2.5 pt 7 "
            "title 'Strassen O(n^{2.81})'\n"
        "unset arrow\nunset label\n"

        // Escala log-log
        "set output '%s'\n"
        "set title 'Estandar vs Strassen - escala log-log\\n"
            "umbral=%d - AMD Ryzen 5 5600H - Ubuntu - g++ -O2 - 15 reps - %.2fs' "
            "font 'Monospace Bold,11'\n"
        "set xlabel 'n (escala log base 2)' font 'Monospace,11'\n"
        "set ylabel 'Tiempo ms (escala log)' font 'Monospace,11'\n"
        "set logscale xy 2\n"
        "set format x '2^{%%L}'\n"
        "set format y '%%.4g'\n"
        "set xtics auto rotate by -45\n"
        "%s"
        "plot '%s' using 1:2 with linespoints lc rgb '#185FA5' lw 2.5 pt 5 "
            "title 'Estandar O(n^3)',\\\n"
        "     '%s' using 1:3 with linespoints lc rgb '#D85A30' lw 2.5 pt 7 "
            "title 'Strassen O(n^{2.81})'\n"
        "unset logscale\nunset arrow\nunset label\n",

        lineal, umbral, totalSeg, marcaN0, csvMedianas, csvMedianas,
        loglog, umbral, totalSeg, marcaN0, csvMedianas, csvMedianas
    );

    FILE* gp = popen("gnuplot", "w");
    if (!gp) throw runtime_error("gnuplot no encontrado. Instalar: sudo apt install gnuplot");
    fputs(script, gp);
    pclose(gp);

    printf("\n  Graficos guardados en:\n    %s\n    %s\n", lineal, loglog);
}

int main() {
    const vector<int> dimensiones = {2, 4, 8, 16, 24, 28, 32, 36, 40, 48, 56, 64, 128, 256, 512, 1024};
    const int REPS   = 15;
    const int UMBRAL = 32;

    fs::create_directories("data");
    fs::create_directories("figures");

    FILE* csvRaw = fopen("data/tiempos.csv", "w");
    if (!csvRaw) {
        fprintf(stderr, "Error: no se pudo crear data/tiempos.csv\n");
        return 1;
    }
    fprintf(csvRaw, "n,algoritmo,tiempo_ms\n");

    vector<ResultadoN> resultados;
    resultados.reserve(dimensiones.size());

    auto t0Total = chrono::high_resolution_clock::now();

    for (int n : dimensiones) {
        printf("  n = %5d  ...  ", n);
        fflush(stdout);
        ResultadoN r = ejecutarBenchmark(n, REPS, UMBRAL, csvRaw);
        resultados.push_back(r);
        printf("Est: %.4f ms  |  Str: %.4f ms\n", r.medianaEstandar, r.medianaStrassen);
    }
    fclose(csvRaw);

    double totalSeg = chrono::duration<double>(
        chrono::high_resolution_clock::now() - t0Total).count();

    imprimirTabla(resultados, UMBRAL);

    int n0 = detectarN0(resultados);
    if (n0 != -1)
        printf("\n  N0 detectado empiricamente: n = %d\n", n0);
    else
        printf("\n  N0 no detectado en este rango.\n"
               "  El overhead de alocacion supera la ganancia teorica en n <= 1024.\n");

    escribirCSVMedianas(resultados, "data/medianas.csv");
    printf("  CSV raw     ->  data/tiempos.csv\n"
           "  CSV resumen ->  data/medianas.csv\n"
           "  Tiempo total: %.2f s\n", totalSeg);

    graficar("data/medianas.csv", "figures", n0, UMBRAL, totalSeg);

    return 0;
}
