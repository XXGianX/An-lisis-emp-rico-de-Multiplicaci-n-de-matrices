# Análisis Empírico: Algoritmo de Strassen vs. Multiplicación Estándar

Este repositorio contiene el desarrollo y análisis del experimento realizado para contrastar la eficiencia teórica del algoritmo de Strassen ($O(n^{2.81})$) frente al algoritmo estándar ($O(n^3)$), buscando determinar el valor crítico $N_0$ a partir del cual Strassen es superior.

## Objetivo del Experimento
Identificar de manera empírica el valor de $N_0$ mediante mediciones controladas de tiempo de ejecución, analizando factores como la jerarquía de memoria y el sobrecosto de la recursividad.

## Especificaciones del Entorno
*   **Lenguaje:** C++17
*   **Compilador:** g++ (versión por definir)
*   **Sistema Operativo:** (por definir)
*   **Hardware:** (por definir)
*   **Herramienta de Graficación:** gnuplot

## Estructura del Proyecto
*   `src/`: Archivos fuente (.cpp, .h).
*   `data/`: Archivos .dat generados con los tiempos de ejecución.
*   `plots/`: Gráficos generados por gnuplot.
*   `Makefile`: Automatización de compilación y ejecución.

## Cómo reproducir
1. Asegúrese de tener instalado `gnuplot`.
2. Ejecute `make` para compilar.
3. Ejecute `./matrix_analysis` para generar datos y gráficas.
