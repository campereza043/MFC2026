#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std;

int main(){

    // parametros
    double rb = 0.002;
    double h = 1.0;
    int N = 50;

    double p0 = 0.01;

    // arreglos
    double t[51];
    double p_exact[51];
    double p_euler[51];
    double p_taylor[51];
    double p_trap[51];

    // inicializacion tiempo y exacta
    for(int i=0; i<=N; i++){
        t[i] = i*h;
        p_exact[i] = 1 - 0.99*exp(-rb*t[i]);
    }

    // condiciones iniciales
    p_euler[0] = p0;
    p_taylor[0] = p0;
    p_trap[0] = p0;

    // =========================
    // EULER
    // =========================

    auto start_euler = chrono::high_resolution_clock::now();

    for(int n=0; n<N; n++)
        p_euler[n+1] = p_euler[n] + rb*(1 - p_euler[n]);

    auto end_euler = chrono::high_resolution_clock::now();
    chrono::duration<double> time_euler = end_euler - start_euler;

    // =========================
    // TAYLOR ORDEN 2
    // =========================

    auto start_taylor = chrono::high_resolution_clock::now();

    for(int n=0; n<N; n++)
        p_taylor[n+1] = p_taylor[n] + 0.001998*(1 - p_taylor[n]);

    auto end_taylor = chrono::high_resolution_clock::now();
    chrono::duration<double> time_taylor = end_taylor - start_taylor;

    // =========================
    // TRAPECIO
    // =========================

    auto start_trap = chrono::high_resolution_clock::now();

    for(int n=0; n<N; n++)
        p_trap[n+1] = (0.999*p_trap[n] + 0.002)/1.001;

    auto end_trap = chrono::high_resolution_clock::now();
    chrono::duration<double> time_trap = end_trap - start_trap;

    // =========================
    // ERRORES ABSOLUTOS EN t=50
    // =========================

    double error_euler = fabs(p_exact[N] - p_euler[N]);
    double error_taylor = fabs(p_exact[N] - p_taylor[N]);
    double error_trap = fabs(p_exact[N] - p_trap[N]);

    // =========================
    // GUARDAR ARCHIVOS .dat
    // =========================

    ofstream exact_file("exacta.dat");
    ofstream euler_file("euler.dat");
    ofstream taylor_file("taylor.dat");
    ofstream trap_file("trapecio.dat");

    for(int i=0; i<=N; i++){
        exact_file << t[i] << " " << p_exact[i] << endl;
        euler_file << t[i] << " " << p_euler[i] << endl;
        taylor_file << t[i] << " " << p_taylor[i] << endl;
        trap_file << t[i] << " " << p_trap[i] << endl;
    }

    exact_file.close();
    euler_file.close();
    taylor_file.close();
    trap_file.close();

    // =========================
    // RESULTADOS EN CONSOLA
    // =========================

    cout << "========== RESULTADOS ==========\n";

    cout << "\n--- Euler ---\n";
    cout << "Valor en t=50: " << p_euler[N] << endl;
    cout << "Error absoluto: " << error_euler << endl;
    cout << "Tiempo ejecucion: " << time_euler.co_