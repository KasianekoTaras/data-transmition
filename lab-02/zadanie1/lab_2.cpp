#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

vector<complex<double>> dft(const vector<double>& x) {
    int N = x.size();
    vector<complex<double>> X(N);
    for (int k = 0; k < N; ++k) {
        complex<double> sum(0, 0);
        for (int n = 0; n < N; ++n) {
            double с1 = (2.0 * M_PI * k * n) / N;
            sum += x[n] * complex<double>(cos(с1), -sin(с1));
        }
        X[k] = sum;
    }
    return X;
}

int main() {
    const double fs = 16000.0;
    const double Tc = 1;
    const int N = static_cast<int>(fs * Tc);
    const double A = 1.0;
    const double f = 1000.0;

    vector<double> t(N), x(N);
    double dt = 1.0 / fs;
    for (int i = 0; i < N; ++i) {
        t[i] = i * dt;
        //x[i] = sin(2.0 * M_PI * 1000.0 * t[i]) + 0.5 * sin(2.0 * M_PI * 500.0 * t[i]);
        x[i] = A * sin(2.0 * M_PI * f * t[i]);
    }

    vector<complex<double>> Xk = dft(x);

    int halfN = N / 2;
    vector<double> fk(halfN), M(halfN);

    for (int k = 0; k < halfN; ++k) {
        fk[k] = k * (fs / N);
        M[k] = sqrt(pow(Xk[k].real(), 2) + pow(Xk[k].imag(), 2));
    }

    plt::subplot(2, 1, 1); plt::plot(t, x); plt::xlim(0.0, 0.03); plt::title("Sygnal w dziedzinie czasu x(t)"); plt::grid(true);
    plt::subplot(2, 1, 2); plt::plot(fk, M); plt::title("Widmo amplitudowe M(k)"); plt::grid(true);

    plt::tight_layout();
    plt::save("zadanie_1.png");
    plt::show();

    return 0;
}