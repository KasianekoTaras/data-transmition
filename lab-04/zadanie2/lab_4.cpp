#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <string>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

void fft(std::vector<std::complex<double>>& a) {
    int n = a.size();
    if (n <= 1) return;

    std::vector<std::complex<double>> a_even(n / 2), a_odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        a_even[i] = a[i * 2];
        a_odd[i] = a[i * 2 + 1];
    }

    fft(a_even);
    fft(a_odd);

    for (int i = 0; i < n / 2; i++) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * i / n) * a_odd[i];
        a[i] = a_even[i] + t;
        a[i + n / 2] = a_even[i] - t;
    }
}

std::vector<double> get_spectrum_db(const std::vector<double>& signal) {
    int n = signal.size();
    std::vector<std::complex<double>> comp_sig(n);
    for (int i = 0; i < n; i++) {
        comp_sig[i] = signal[i];
    }

    fft(comp_sig);

    std::vector<double> mag_db(n / 2);
    for (int i = 0; i < n / 2; i++) {
        double mag = std::abs(comp_sig[i]) / (n / 2.0);
        mag_db[i] = 20 * std::log10(mag + 1e-12);
    }
    return mag_db;
}

int main() {
    double fm = 2.0;
    double fn = 50.0;
    double fs = 1024.0;
    double duration = 1.0;
    double dt = 1.0 / fs;
    int num_samples = duration * fs;

    std::vector<double> t(num_samples);
    std::vector<double> m(num_samples);
    std::vector<double> freqs(num_samples / 2);

    for (int i = 0; i < num_samples; ++i) {
        t[i] = i * dt;
        m[i] = std::sin(2 * M_PI * fm * t[i]);
    }

    for (int i = 0; i < num_samples / 2; ++i) {
        freqs[i] = i * fs / num_samples;
    }
    int plot_idx = 1;

    std::vector<double> kA_vals = { 0.5, 5.0, 25.0 };
    for (double kA : kA_vals) {
        std::vector<double> zA(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            zA[i] = (kA * m[i] + 1.0) * std::cos(2 * M_PI * fn * t[i]);
        }
        std::vector<double> MA = get_spectrum_db(zA);

        plt::subplot(3, 3, plot_idx++);
        plt::plot(freqs, MA);
        plt::xlim(0, 100);
        plt::ylim(-80, 10);
        plt::title("Widmo AM: kA = " + std::to_string(kA));
    }

    std::vector<double> kP_vals = { 0.5, 2.0, 8.0 };
    for (double kP : kP_vals) {
        std::vector<double> zP(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            zP[i] = std::cos(2 * M_PI * fn * t[i] + kP * m[i]);
        }
        std::vector<double> MP = get_spectrum_db(zP);

        plt::subplot(3, 3, plot_idx++);
        plt::plot(freqs, MP);
        plt::xlim(0, 100);
        plt::ylim(-80, 10);
        plt::title("Widmo PM: kP = " + std::to_string(kP));
    }

    std::vector<double> kF_vals = { 0.5, 2.0, 8.0 };
    for (double kF : kF_vals) {
        std::vector<double> zF(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            zF[i] = std::cos(2 * M_PI * fn * t[i] + (kF / fm) * m[i]);
        }
        std::vector<double> MF = get_spectrum_db(zF);

        plt::subplot(3, 3, plot_idx++);
        plt::plot(freqs, MF);
        plt::xlim(0, 100);
        plt::ylim(-80, 10);
        plt::title("Widmo FM: kF = " + std::to_string(kF));
    }

    plt::tight_layout();
    plt::save("zadanie_1.png");
    plt::show();

    return 0;
}