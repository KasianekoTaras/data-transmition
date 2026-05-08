#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <fftw3.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

vector<int> stringToBits(const string& text) {
    vector<int> bits;
    for (char c : text) {
        if (c >= 32 && c <= 127) {
            for (int i = 6; i >= 0; --i) {
                bits.push_back((c >> i) & 1);
            }
        }
    }
    return bits;
}

vector<double> getSpectrumDecibels(const vector<double>& signal) {
    int N = signal.size();
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    double* in = (double*)fftw_malloc(sizeof(double) * N);
    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) in[i] = signal[i];
    fftw_execute(p);

    vector<double> spectrum(N / 2 + 1);
    for (int i = 0; i < N / 2 + 1; ++i) {
        double mag = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        mag = max(mag, 1e-12);
        spectrum[i] = max(0.0, 20.0 * log10(mag));
    }
    fftw_destroy_plan(p); fftw_free(in); fftw_free(out);
    return spectrum;
}

void calculateBandwidthDB(const vector<double>& f, const vector<double>& spectrum, double drop_dB, double& f_min, double& f_max, double& B) {
    double max_val = *max_element(spectrum.begin(), spectrum.end());
    double threshold = max_val - drop_dB;

    f_min = f.back();
    f_max = f.front();

    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (spectrum[i] >= threshold) {
            if (f[i] < f_min) f_min = f[i];
            if (f[i] > f_max) f_max = f[i];
        }
    }
    B = f_max - f_min;
}

string formatDouble(double value) {
    ostringstream out;
    out << fixed << setprecision(1) << value;
    return out.str();
}

int main() {
    string text = "DATA";
    vector<int> b = stringToBits(text);
    int B_bits = b.size();

    double Tb = 0.1;
    double Tc = B_bits * Tb;
    int W = 2;
    double fn = W / Tb;
    double fn1 = (W + 1) / Tb;
    double fn2 = (W + 2) / Tb;

    double fs = 4.0 * fn2 * 10;
    double dt = 1.0 / fs;
    int N = ceil(Tc / dt);

    vector<double> t(N), zA(N), zP(N), zF(N);
    for (int i = 0; i < N; ++i) {
        t[i] = i * dt;
        int bit_idx = t[i] / Tb;
        if (bit_idx >= B_bits) bit_idx = B_bits - 1;

        if (b[bit_idx] == 0) {
            zA[i] = 0.5 * sin(2.0 * M_PI * fn * t[i]);
            zP[i] = sin(2.0 * M_PI * fn * t[i]);
            zF[i] = sin(2.0 * M_PI * fn1 * t[i]);
        }
        else {
            zA[i] = 1.0 * sin(2.0 * M_PI * fn * t[i]);
            zP[i] = sin(2.0 * M_PI * fn * t[i] + M_PI);
            zF[i] = sin(2.0 * M_PI * fn2 * t[i]);
        }
    }

    vector<double> MA = getSpectrumDecibels(zA);
    vector<double> MP = getSpectrumDecibels(zP);
    vector<double> MF = getSpectrumDecibels(zF);

    vector<double> f(MA.size());
    for (size_t i = 0; i < f.size(); ++i) f[i] = i * fs / N;

    ofstream outfile("zadanie_1.txt");
    double fmin, fmax, bw;
    vector<string> names = { "ASK", "PSK", "FSK" };
    vector<vector<double>> spectra = { MA, MP, MF };
    vector<double> drops = { 3.0, 6.0, 10.0 };
    vector<string> colors = { "r", "g", "b" };

    outfile << "\n\n";
    for (size_t s = 0; s < 3; ++s) {
        outfile << names[s] << "\n";
        for (double drop : drops) {
            calculateBandwidthDB(f, spectra[s], drop, fmin, fmax, bw);
            outfile << "B_" << drop << "dB = " << bw << " Hz\n";
        }
        outfile << "\n";
    }

    outfile << "\n";
    outfile << "1. Wraz ze zwiekszeniem analizowanego spadku poziomu odniesienia (z 3 dB do 10 dB) wyznaczana szerokosc pasma rosnie. Wynika to z faktu uwzgledniania coraz nizszych, szerszych partii widma.\n";
    outfile << "2. Rozne typy modulacji charakteryzuja sie odmienna szerokoscia pasma. Modulacja FSK zazwyczaj zajmuje szersze pasmo w porownaniu do ASK i PSK ze wzgledu na obecnosc dwoch roznych czestotliwosci nosnych w sygnale.\n";
    outfile.close();

    for (int i = 0; i < 3; ++i) {
        plt::subplot(3, 1, i + 1);
        plt::stem(f, spectra[i], "k");
        plt::plot(f, spectra[i], "k.");

        double max_val = *max_element(spectra[i].begin(), spectra[i].end());

        plt::plot({ 0.0, 100.0 }, { max_val, max_val }, "k--");

        for (size_t d = 0; d < drops.size(); ++d) {
            double drop = drops[d];
            string c = colors[d];
            double threshold = max_val - drop;

            calculateBandwidthDB(f, spectra[i], drop, fmin, fmax, bw);

            plt::plot({ 0.0, 100.0 }, { threshold, threshold }, c + "--");
            plt::plot({ fmin, fmin }, { 0.0, threshold }, c + "-");
            plt::plot({ fmax, fmax }, { 0.0, threshold }, c + "-");
        }

        plt::title("Widmo amplitudowe " + names[i] + " [dB]");
        if (i == 2) plt::xlabel("Czestotliwosc [Hz]");
        plt::ylabel("Amplituda [dB]");
        plt::xlim(0.0, 100.0);
        plt::ylim(0.0, max_val + 15.0);
        plt::grid(true);
    }

    plt::save("zadanie_1.png");
    plt::show();

    return 0;
}