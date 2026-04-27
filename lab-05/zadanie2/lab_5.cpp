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
#include <sstream>
#include <fftw3.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

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
    fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);
    return mag_lin;
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

        std::vector<double> M_lin = get_spectrum_lin(signal);

        double E_total = 0;
        for (double val : M_lin) E_total += val * val;

        double alpha = 0.0, r_alpha = 0.0, E_alpha = 0.0;

        std::ostringstream log_stream;
        log_stream << std::fixed << std::setprecision(2);

        while (true) {
            alpha += 1.0;
            int ia = std::round(((fn - alpha) * N) / fs);
            int ib = std::round(((fn + alpha) * N) / fs);
            if (ia < 0) ia = 0;
            if (ib >= (int)M_lin.size()) ib = M_lin.size() - 1;

            E_alpha = 0;
            for (int k = ia; k <= ib; ++k) E_alpha += M_lin[k] * M_lin[k];

            r_alpha = (E_alpha / E_total) * 100.0;

            log_stream << "    -> alpha = " << std::setw(4) << alpha << " Hz | E_alpha = "
                << std::setw(6) << E_alpha << " | r_alpha = "
                << std::setw(6) << r_alpha << "%\n";

            if (r_alpha > 80.0 || ib == (int)M_lin.size() - 1) break;
        }

        double B = 2.0 * alpha;

        outfile << "Sygnal: " << names[j] << " (k = " << std::fixed << std::setprecision(2) << k_vals[j] << ")\n"
            << "  Energia calkowita (E) = " << E_total << "\n"
            << "  Energia w oknie (E_alpha) = " << E_alpha << " (r_alpha = " << r_alpha << "%)\n"
            << "  Szerokosc pasma B = " << B << " Hz (alpha = " << alpha << " Hz)\n"
            << "  Przebieg obliczen:\n"
            << log_stream.str();

        plt::subplot(3, 3, j + 1);
        std::vector<double> sx, sy, px, py;
        for (size_t i = 0; i < freqs.size(); ++i) {
            sx.push_back(freqs[i]); sy.push_back(0.0);
            sx.push_back(freqs[i]); sy.push_back(M_lin[i]);
            sx.push_back(std::numeric_limits<double>::quiet_NaN());
            sy.push_back(std::numeric_limits<double>::quiet_NaN());

            if (M_lin[i] > 0.001) {
                px.push_back(freqs[i]);
                py.push_back(M_lin[i]);
            }
        }
        plt::plot(sx, sy, "C0-");
        plt::plot(px, py, "C0o");

        double max_lin = *std::max_element(M_lin.begin(), M_lin.end());
        plt::plot({ fn - alpha, fn - alpha }, { 0.0, max_lin }, "g--");
        plt::plot({ fn + alpha, fn + alpha }, { 0.0, max_lin }, "g--");

        plt::title(names[j] + " k=" + std::to_string(k_vals[j]).substr(0, 4));
        plt::xlim(0.0, 100.0);
    }

    outfile << "\nWNIOSKI I OBSERWACJE:\n"
        << "1. Szerokosc pasma wyznaczona metoda energetyczna skupia sie na obszarze zawierajacym 80% mocy sygnalu.\n"
        << "2. Metoda ta pozwala na bardziej obiektywne szacowanie pasma dla sygnalow o skomplikowanym widmie (PM, FM).\n"
        << "3. Zwiekszenie wspolczynnika k w modulacjach katowych prowadzi do wyraznego wzrostu szerokosci pasma B.\n";

    plt::tight_layout();
    plt::save("zadanie_2.png");
    plt::show();
    outfile.close();
    return 0;
}