    /*
     * * 1. Czestotliwosc Nyquista (granica poprawnego odtwarzania bez aliasingu):
     * f_max = fs / 2 = 512.0 / 2 = 256.0 Hz
     * * 2. Maksymalna teoretyczna liczba prazkow (harmonicznych) mieszczacych sie w pasmie (0, 256 Hz):
     * * - Dla sygnalu x(t) [piloksztaltny - zawiera wszystkie wielokrotnosci f]:
     * Wzor: k * f < f_max
     * k * 2.0 < 256.0
     * k < 128
     * Wniosek: Maksymalnie 127 prazkow.
     * * - Dla sygnalu y(t) i z(t) [trojkatny/prostokatny - zawiera tylko nieparzyste wielokrotnosci f]:
     * Wzor: (2k - 1) * f < f_max
     * (2k - 1) * 2.0 < 256.0
     * 2k - 1 < 128
     * 2k < 129
     * k < 64.5
     * Wniosek: Maksymalnie 64 prazki.
     * * - Sygnal x(t) uzywa H = 5 (5 kolejnych harmonicznych):
     * Czestotliwosci (k * f): 1*2=2, 2*2=4, 3*2=6, 4*2=8, 5*2=10 Hz.
     * Na wykresie (d) bedzie widoczne dokladnie 5 prazkow (do 10 Hz).
     * - Sygnaly y(t) i z(t) uzywaja H = 5 (5 nieparzystych harmonicznych):
     * Czestotliwosci ((2k-1) * f): 1*2=2, 3*2=6, 5*2=10, 7*2=14, 9*2=18 Hz.
     * Na wykresach (e) i (f) bedzie widocznych 5 prazkow (do 18 Hz).
     */

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    double f = 2.0;
    double fs = 512.0;
    int H = 5;
    int N = 1024.0;
    double dt = 1.0 / fs;

    std::vector<double> t_vec(N), x(N, 0.0), y(N, 0.0), z(N, 0.0);

    fftw_complex* in_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* out_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan_x = fftw_plan_dft_1d(N, in_x, out_x, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_complex* in_y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* out_y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan_y = fftw_plan_dft_1d(N, in_y, out_y, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_complex* in_z = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* out_z = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan_z = fftw_plan_dft_1d(N, in_z, out_z, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) {
        double t = i * dt;
        t_vec[i] = t;

        for (int k = 1; k <= H; ++k) {
            x[i] += pow(-1, k + 1) * sin(2 * M_PI * k * f * t) / k;
        }
        x[i] *= (2.0 / M_PI);

        for (int k = 1; k <= H; ++k) {
            int odd_k = 2 * k - 1;
            y[i] += pow(-1, k - 1) * sin(2 * M_PI * odd_k * f * t) / (odd_k * odd_k);
            z[i] += sin(2 * M_PI * odd_k * f * t) / odd_k;
        }
        y[i] *= (8.0 / (M_PI * M_PI));
        z[i] *= (4.0 / M_PI);

        in_x[i][0] = x[i]; in_x[i][1] = 0.0;
        in_y[i][0] = y[i]; in_y[i][1] = 0.0;
        in_z[i][0] = z[i]; in_z[i][1] = 0.0;
    }

    fftw_execute(plan_x);
    fftw_execute(plan_y);
    fftw_execute(plan_z);

    std::vector<double> fk, Mx, My, Mz;
    std::vector<double> fk_x, val_x, fk_y, val_y, fk_z, val_z;

    for (int i = 0; i < N / 2; ++i) {
        double current_fk = i * (fs / N);
        double current_Mx = 2.0 * sqrt(out_x[i][0] * out_x[i][0] + out_x[i][1] * out_x[i][1]) / N;
        double current_My = 2.0 * sqrt(out_y[i][0] * out_y[i][0] + out_y[i][1] * out_y[i][1]) / N;
        double current_Mz = 2.0 * sqrt(out_z[i][0] * out_z[i][0] + out_z[i][1] * out_z[i][1]) / N;

        fk.push_back(current_fk);
        Mx.push_back(current_Mx);
        My.push_back(current_My);
        Mz.push_back(current_Mz);

        if (current_Mx > 1e-3) { fk_x.push_back(current_fk); val_x.push_back(current_Mx); }
        if (current_My > 1e-3) { fk_y.push_back(current_fk); val_y.push_back(current_My); }
        if (current_Mz > 1e-3) { fk_z.push_back(current_fk); val_z.push_back(current_Mz); }
    }

    fftw_destroy_plan(plan_x); fftw_free(in_x); fftw_free(out_x);
    fftw_destroy_plan(plan_y); fftw_free(in_y); fftw_free(out_y);
    fftw_destroy_plan(plan_z); fftw_free(in_z); fftw_free(out_z);


    plt::subplot(3, 2, 1); plt::plot(t_vec, x); plt::xlim(0.0, 1.5);
    plt::subplot(3, 2, 2); plt::stem(fk_x, val_x, "k"); plt::plot({ 0.0, 20.0 }, { 0.0, 0.0 }, "k"); plt::xlim(0.0, 20.0); plt::xlabel("(d)");

    plt::subplot(3, 2, 3); plt::plot(t_vec, y); plt::xlim(0.0, 1.5);
    plt::subplot(3, 2, 4); plt::stem(fk_y, val_y, "k"); plt::plot({ 0.0, 20.0 }, { 0.0, 0.0 }, "k"); plt::xlim(0.0, 20.0); plt::xlabel("(e)");

    plt::subplot(3, 2, 5); plt::plot(t_vec, z); plt::xlim(0.0, 1.5);
    plt::subplot(3, 2, 6); plt::stem(fk_z, val_z, "k"); plt::plot({ 0.0, 20.0 }, { 0.0, 0.0 }, "k"); plt::xlim(0.0, 20.0); plt::xlabel("(f)");

    plt::tight_layout();
    plt::save("zadanie_1.png");
    plt::show();

    return 0;
}