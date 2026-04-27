#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <fftw3.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::vector<int> stringToBits(const std::string& text) {
    std::vector<int> bits;
    for (char c : text) {
        if (c >= 32 && c <= 127) {
            for (int i = 6; i >= 0; --i) {
                bits.push_back((c >> i) & 1);
            }
        }
    }
    return bits;
}

std::vector<double> getSpectrumDecibels(const std::vector<double>& signal) {
    int N = signal.size();
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    double* in = (double*)fftw_malloc(sizeof(double) * N);

    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) {
        in[i] = signal[i];
    }

    fftw_execute(p);

    std::vector<double> spectrum(N / 2 + 1);
    for (int i = 0; i < N / 2 + 1; ++i) {
        double mag = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / (N / 2.0);
        mag = std::max(mag, 1e-12);
        spectrum[i] = 20.0 * log10(mag);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return spectrum;
}

int main() {
    std::string text = "DATA";
    std::vector<int> b = stringToBits(text);
    int B = b.size();

    double Tb = 0.1;
    double Tc = B * Tb;

    int W = 2;
    double fn = W / Tb;
    double fn1 = (W + 1) / Tb;
    double fn2 = (W + 2) / Tb;

    double fs = 4.0 * fn2 * 10;
    double dt = 1.0 / fs;
    int N = std::ceil(Tc / dt);

    double A1 = 0.5;
    double A2 = 1.0;

    std::vector<double> t(N);
    std::vector<double> zA(N), zP(N), zF(N);

    for (int i = 0; i < N; ++i) {
        t[i] = i * dt;
        int bit_idx = t[i] / Tb;
        if (bit_idx >= B) bit_idx = B - 1;
        int bit = b[bit_idx];

        if (bit == 0) {
            zA[i] = A1 * sin(2.0 * M_PI * fn * t[i]);
            zP[i] = sin(2.0 * M_PI * fn * t[i]);
            zF[i] = sin(2.0 * M_PI * fn1 * t[i]);
        }
        else {
            zA[i] = A2 * sin(2.0 * M_PI * fn * t[i]);
            zP[i] = sin(2.0 * M_PI * fn * t[i] + M_PI);
            zF[i] = sin(2.0 * M_PI * fn2 * t[i]);
        }
    }

    std::vector<double> MA = getSpectrumDecibels(zA);
    std::vector<double> MP = getSpectrumDecibels(zP);
    std::vector<double> MF = getSpectrumDecibels(zF);

    std::vector<double> f(MA.size());
    for (size_t i = 0; i < f.size(); ++i) {
        f[i] = i * fs / N;
    }

    plt::subplot(3, 1, 1);
    plt::plot(f, MA);
    plt::title("Widmo amplitudowe ASK [dB]");
    plt::ylabel("Amplituda [dB]");
    plt::xlim(0.0, 100.0);
    plt::ylim(-100.0, 0.0);
    plt::grid(true);

    plt::subplot(3, 1, 2);
    plt::plot(f, MF);
    plt::title("Widmo amplitudowe FSK [dB]");
    plt::ylabel("Amplituda [dB]");
    plt::xlim(0.0, 100.0);
    plt::ylim(-100.0, 0.0);
    plt::grid(true);

    plt::subplot(3, 1, 3);
    plt::plot(f, MP);
    plt::title("Widmo amplitudowe PSK [dB]");
    plt::xlabel("Czestotliwosc [Hz]");
    plt::ylabel("Amplituda [dB]");
    plt::xlim(0.0, 100.0);
    plt::ylim(-100.0, 0.0);
    plt::grid(true);

    plt::save("widma_db.png");
    plt::show();

    return 0;
}