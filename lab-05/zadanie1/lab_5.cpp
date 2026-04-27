#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <limits>
#include <fftw3.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

struct BandwidthResult {
    double B;
    double f_min;
    double f_max;
};

BandwidthResult calculate_bandwidth(const std::vector<double>& freqs, const std::vector<double>& mag_db, double drop) {
    auto max_it = std::max_element(mag_db.begin(), mag_db.end());
    double max_val = *max_it;
    double threshold = max_val - drop;

    int min_idx = 0;
    int max_idx_cross = 0;

    for (int i = 0; i < mag_db.size(); ++i) {
        if (mag_db[i] >= threshold) {
            min_idx = i;
            break;
        }
    }
    for (int i = mag_db.size() - 1; i >= 0; --i) {
        if (mag_db[i] >= threshold) {
            max_idx_cross = i;
            break;
        }
    }
    return { freqs[max_idx_cross] - freqs[min_idx], freqs[min_idx], freqs[max_idx_cross] };
}

std::vector<double> get_spectrum_lin(const std::vector<double>& signal) {
    int n = signal.size();
    double* in = (double*)fftw_malloc(sizeof(double) * n);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1));
    for (int i = 0; i < n; ++i) in[i] = signal[i];

    fftw_plan plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    std::vector<double> mag_lin(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        mag_lin[i] = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / (n / 2.0);
    }
    fftw_destroy_plan(plan); 
    fftw_free(in); 
    fftw_free(out);
    return mag_lin;
}

std::vector<double> get_amplitude_spectrum_db(const std::vector<double>& signal) {
    std::vector<double> mag_lin = get_spectrum_lin(signal);
    std::vector<double> mag_db(mag_lin.size());
    for (size_t i = 0; i < mag_lin.size(); i++) mag_db[i] = 20 * std::log10(mag_lin[i] + 1e-12);
    return mag_db;
}

int main() {
    double fm = 2.0, fn = 50.0, fs = 1024.0, duration = 1.0;
    double dt = 1.0 / fs;
    int N = duration * fs;

    std::vector<double> t(N), m(N), freqs(N / 2);
    for (int i = 0; i < N; ++i) {
        t[i] = i * dt;
        m[i] = std::sin(2 * M_PI * fm * t[i]);
    }
    for (int i = 0; i < N / 2; ++i) freqs[i] = i * fs / N;

    std::vector<std::string> names = { "AM", "AM", "AM", "PM", "PM", "PM", "FM", "FM", "FM" };
    std::vector<double> k_vals = { 0.5, 5.0, 25.0, 0.5, 2.0, 8.0, 0.5, 2.0, 8.0 };

    std::ofstream outfile("wyniki.txt");
    outfile << "Parametry: fm=2Hz, fn=50Hz, fs=1024Hz\n\n";

    for (int j = 0; j < 9; ++j) {
        std::vector<double> signal(N);
        if (names[j] == "AM") for (int i = 0; i < N; i++) signal[i] = (k_vals[j] * m[i] + 1.0) * std::cos(2 * M_PI * fn * t[i]);
        else if (names[j] == "PM") for (int i = 0; i < N; i++) signal[i] = std::cos(2 * M_PI * fn * t[i] + k_vals[j] * m[i]);
        else if (names[j] == "FM") for (int i = 0; i < N; i++) signal[i] = std::cos(2 * M_PI * fn * t[i] + (k_vals[j] / fm) * m[i]);

        std::vector<double> M_db = get_amplitude_spectrum_db(signal);
        BandwidthResult b3 = calculate_bandwidth(freqs, M_db, 3.0);
        BandwidthResult b6 = calculate_bandwidth(freqs, M_db, 6.0);
        BandwidthResult b10 = calculate_bandwidth(freqs, M_db, 10.0);

        outfile << "Sygnal: " << names[j] << " (k = " << k_vals[j] << ")\n" << std::fixed << std::setprecision(2)
            << "  B_3dB  = " << b3.B << " Hz (" << b3.f_min << "-" << b3.f_max << " Hz)\n"
            << "  B_6dB  = " << b6.B << " Hz (" << b6.f_min << "-" << b6.f_max << " Hz)\n"
            << "  B_10dB = " << b10.B << " Hz (" << b10.f_min << "-" << b10.f_max << " Hz)\n";

        plt::subplot(3, 3, j + 1);
        double max_v = *std::max_element(M_db.begin(), M_db.end()), bot = max_v - 50.0;
        std::vector<double> sx, sy;
        for (size_t i = 0; i < freqs.size(); ++i) {
            sx.push_back(freqs[i]); sy.push_back(bot); sx.push_back(freqs[i]); sy.push_back(M_db[i]);
            sx.push_back(std::numeric_limits<double>::quiet_NaN()); sy.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        plt::plot(sx, sy, "C0-"); 
        plt::plot(freqs, M_db, "C0o");
        plt::plot({ 0.0, 100.0 }, { max_v, max_v }, "r--");
        plt::plot({ 0.0, 100.0 }, { max_v - 3.0, max_v - 3.0 }, "g:");
        plt::plot({ 0.0, 100.0 }, { max_v - 6.0, max_v - 6.0 }, "y:");
        plt::plot({ 0.0, 100.0 }, { max_v - 10.0, max_v - 10.0 }, "r--");
        plt::title(names[j] + " k=" + std::to_string(k_vals[j]).substr(0, 4));
        plt::xlim(0.0, 100.0); 
        plt::ylim(bot, max_v + 5.0);
    }

    outfile << "\nWNIOSKI I OBSERWACJE:\n"
        << "1. Szerokosc pasma sygnalu modulowanego zalezy od wspolczynnika modulacji k.\n"
        << "2. Dla modulacji AM szerokosc pasma teoretycznie wynosi 2*fm (4Hz), co widac przy spadku 10dB.\n"
        << "3. W modulacjach PM i FM zwiekszenie k powoduje znaczne rozszerzenie widma (powstawanie bocznych prazkow).\n"
        << "4. Szacowanie pasma zalezy od przyjetego kryterium (3dB, 6dB czy 10dB) - im wiekszy spadek, tym wieksza zmierzona szerokosc.\n";

    plt::tight_layout();
    plt::save("zadanie_1.png");
    plt::show();
    outfile.close();
    return 0;
}