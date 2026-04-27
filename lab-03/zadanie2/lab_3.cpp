#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    double f1 = 5.0;
    double f2 = 12.0;
    double alpha = 2.0;
    double beta = 3.0;
    double fs = 1024.0;
    int N = 1024;
    double dt = 1.0 / fs;

    std::vector<double> t_vec(N), x(N), y(N), z(N);

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

        x[i] = 0.5 * sin(2 * M_PI * f1 * t);
        y[i] = sin(2 * M_PI * f2 * t) + 0.7 * sin(2 * M_PI * f1 * t);
        z[i] = alpha * x[i] + beta * y[i];

        in_x[i][0] = x[i]; in_x[i][1] = 0.0;
        in_y[i][0] = y[i]; in_y[i][1] = 0.0;
        in_z[i][0] = z[i]; in_z[i][1] = 0.0;
    }

    fftw_execute(plan_x);
    fftw_execute(plan_y);
    fftw_execute(plan_z);

    std::vector<double> Mx(N / 2), My(N / 2), Mz(N / 2), Mz_hat(N / 2), fk(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        fk[i] = i * (fs / N);
        Mx[i] = 2.0 * sqrt(out_x[i][0] * out_x[i][0] + out_x[i][1] * out_x[i][1]) / N;
        My[i] = 2.0 * sqrt(out_y[i][0] * out_y[i][0] + out_y[i][1] * out_y[i][1]) / N;
        Mz[i] = 2.0 * sqrt(out_z[i][0] * out_z[i][0] + out_z[i][1] * out_z[i][1]) / N;
        Mz_hat[i] = alpha * Mx[i] + beta * My[i];
    }

    fftw_destroy_plan(plan_x); fftw_free(in_x); fftw_free(out_x);
    fftw_destroy_plan(plan_y); fftw_free(in_y); fftw_free(out_y);
    fftw_destroy_plan(plan_z); fftw_free(in_z); fftw_free(out_z);



    plt::subplot(4, 2, 1); plt::plot(t_vec, x); plt::xlim(0.0, 1.0); plt::title("Sygnal x(t)"); plt::grid(true);
    plt::subplot(4, 2, 2); plt::plot(fk, Mx, "bo"); plt::xlim(0.0, 20.0); plt::title("Widmo Mx(k)"); plt::grid(true);

    plt::subplot(4, 2, 3); plt::plot(t_vec, y); plt::xlim(0.0, 1.0); plt::title("Sygnal y(t)"); plt::grid(true);
    plt::subplot(4, 2, 4); plt::plot(fk, My, "bo"); plt::xlim(0.0, 20.0); plt::title("Widmo My(k)"); plt::grid(true);

    plt::subplot(4, 2, 5); plt::plot(t_vec, z); plt::xlim(0.0, 1.0); plt::title("Sygnal z(t)"); plt::grid(true);
    plt::subplot(4, 2, 6); plt::plot(fk, Mz, "bo"); plt::xlim(0.0, 20.0); plt::title("Widmo Mz(k)"); plt::grid(true);

    plt::subplot(4, 2, 8); plt::plot(fk, Mz_hat, "bo"); plt::xlim(0.0, 20.0); plt::title("Widmo Mz_hat"); plt::grid(true);

    plt::tight_layout();
    plt::save("zadanie_2.png");
    plt::show();

    return 0;
}