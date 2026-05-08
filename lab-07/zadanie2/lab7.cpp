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

vector<double> getSpectrumLinear(const vector<double>& signal) {
    int N = signal.size();
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    double* in = (double*)fftw_malloc(sizeof(double) * N);
    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) in[i] = signal[i];
    fftw_execute(p);

    vector<double> spectrum(N / 2 + 1);
    for (int i = 0; i < N / 2 + 1; ++i) {
        spectrum[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / (N / 2.0);
    }
    fftw_destroy_plan(p); fftw_free(in); fftw_free(out);
    return spectrum;
}

void calculateBandwidthEnergy(const vector<double>& X, double fn, double fs, int N_total, double& alpha, double& B_energy, double& r_alpha) {
    double E_total = 0.0;
    for (double mag : X) {
        E_total += mag * mag;
    }

    alpha = 1.0;
    r_alpha = 0.0;

    while (alpha < fs / 2.0) {
        double fa = fn - alpha;
        double fb = fn + alpha;

        int ia = max(0, (int)round((fa * N_total) / fs));
        int ib = min((int)X.size() - 1, (int)round((fb * N_total) / fs));

        double E_alpha = 0.0;
        for (int k = ia; k <= ib; ++k) {
            E_alpha += X[k] * X[k];
        }

        r_alpha = (E_alpha / E_total) * 100.0;

        if (r_alpha > 80.0) {
            B_energy = fb - fa;
            break;
        }
        alpha += 0.5;
    }
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

    vector<double> XA = getSpectrumLinear(zA);
    vector<double> XP = getSpectrumLinear(zP);
    vector<double> XF = getSpectrumLinear(zF);

    vector<double> f(XA.size());
    for (size_t i = 0; i < f.size(); ++i) f[i] = i * fs / N;

    ofstream outfile("zadanie_2.txt");
    double alpha, bw, r_alpha;
    vector<string> names = { "ASK", "PSK", "FSK" };
    vector<vector<double>> spectra = { XA, XP, XF };

    outfile << "\n\n";
    for (size_t s = 0; s < 3; ++s) {
        double current_fn = fn;
        calculateBandwidthEnergy(spectra[s], current_fn, fs, N, alpha, bw, r_alpha);

        outfile << names[s] << "\n";
        outfile << "B = " << bw << " Hz\n";
        outfile << "r_alpha = " << r_alpha << " %\n\n";
    }

    outfile << "\n";
    outfile << "1. Metoda energetyczna pozwala na rzetelne oszacowanie niezbednego pasma transmisyjnego poprzez analize procentowego udzialu calkowitej energii sygnalu (w tym przypadku > 80%).\n";
    outfile << "2. Chociaz widma sygnalow zmodulowanych cyfrowo sa teoretycznie nieograniczone, znakomita wiekszosc energii uzytecznej koncentruje sie w stosunkowo waskim otoczeniu czestotliwosci nosnej, co potwiedza wyznaczone okno (2*alpha).\n";
    outfile.close();

    for (int i = 0; i < 3; ++i) {
        plt::subplot(3, 1, i + 1);
        plt::stem(f, spectra[i], "k");
        plt::plot(f, spectra[i], "k.");

        double max_val = *max_element(spectra[i].begin(), spectra[i].end());
        double current_fn = fn;

        calculateBandwidthEnergy(spectra[i], current_fn, fs, N, alpha, bw, r_alpha);

        double fa = current_fn - alpha;
        double fb = current_fn + alpha;
        double line_top = max_val * 1.1;

        plt::plot({ fa, fa }, { 0.0, line_top }, "g--");
        plt::plot({ fb, fb }, { 0.0, line_top }, "g--");

        plt::title("Widmo liniowe " + names[i]);
        if (i == 2) plt::xlabel("Czestotliwosc [Hz]");
        plt::ylabel("Amplituda");
        plt::xlim(0.0, 100.0);
        plt::ylim(0.0, max_val * 1.2);
        plt::grid(true);
    }

    plt::save("zadanie_2.png");
    plt::show();

    return 0;
}