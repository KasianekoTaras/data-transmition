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
            double angle = 2.0 * M_PI * k * n / N;
            sum += x[n] * complex<double>(cos(angle), -sin(angle));
        }
        X[k] = sum;
    }
    return X;
}

int main() {
    const double fs1 = 16000.0;
    const double Tc1 = 0.05;
    const int N1 = static_cast<int>(fs1 * Tc1);
    const double f_sig1 = 1000.0;

    vector<double> t1(N1), x1(N1);
    double dt1 = 1.0 / fs1;
    for (int i = 0; i < N1; ++i) {
        t1[i] = i * dt1;
        x1[i] = sin(2.0 * M_PI * f_sig1 * t1[i]);
    }

    vector<complex<double>> Xk1 = dft(x1);
    int halfN1 = N1 / 2;
    vector<double> fk1(halfN1), M_dB1(halfN1);

    for (int k = 0; k < halfN1; ++k) {
        fk1[k] = k * (fs1 / N1);
        double Mk = sqrt(pow(Xk1[k].real(), 2) + pow(Xk1[k].imag(), 2));
        double mk_safe = (Mk < 1e-10) ? 1e-10 : Mk;
        M_dB1[k] = 10.0 * log10(mk_safe);
    }

    const double fs2 = 2000.0;
    const double Tc2 = 1.0;
    const int N2 = static_cast<int>(fs2 * Tc2);

    const double f1 = 10.0;
    const double f2 = (fs2 / 2.0) - f1;
    const double f3 = f1 / 2.0;

    vector<double> t2(N2), x2(N2);
    double dt2 = 1.0 / fs2;
    for (int i = 0; i < N2; ++i) {
        t2[i] = i * dt2;
        x2[i] = sin(2.0 * M_PI * f1 * t2[i]) + sin(2.0 * M_PI * f2 * t2[i]) + sin(2.0 * M_PI * f3 * t2[i]);
    }

    vector<complex<double>> Xk2 = dft(x2);
    int halfN2 = N2 / 2;
    vector<double> fk2(halfN2), M_dB2(halfN2);

    for (int k = 0; k < halfN2; ++k) {
        fk2[k] = k * (fs2 / N2);
        double Mk = sqrt(pow(Xk2[k].real(), 2) + pow(Xk2[k].imag(), 2));
        double mk_safe = (Mk < 1e-10) ? 1e-10 : Mk;
        M_dB2[k] = 10.0 * log10(mk_safe);
    }

    plt::subplot(3, 1, 1); plt::plot(fk1, M_dB1); plt::title("Cz. 2 pkt 1: Widmo w skali decybelowej M'(k)"); plt::xlabel("f [Hz]"); plt::ylabel("M'(k) [dB]"); plt::grid(true);
    plt::subplot(3, 1, 2); plt::plot(fk2, M_dB2); plt::title("Cz. 2 pkt 2: 3 tony - Liniowa skala czestotliwosci"); plt::xlabel("f [Hz]"); plt::ylabel("M'(k) [dB]"); plt::grid(true);
    plt::subplot(3, 1, 3); plt::semilogx(fk2, M_dB2); plt::title("Cz. 2 pkt 2: 3 tony - Logarytmiczna skala czestotliwosci"); plt::xlabel("f [Hz]"); plt::ylabel("M'(k) [dB]"); plt::grid(true);

    plt::tight_layout();
    plt::save("zadanie_2.png");
    plt::show();

    return 0;
}