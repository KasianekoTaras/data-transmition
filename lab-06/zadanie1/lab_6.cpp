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

int main() {
    std::string text = "DATA";
    std::vector<int> b = stringToBits(text);
    int B = b.size();

    double Tb = 0.1;
    double Tc = B * Tb;

    int W = 15;
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

    int plot_N = std::min(N, static_cast<int>(3 * Tb / dt));
    std::vector<double> t_plot(t.begin(), t.begin() + plot_N);
    std::vector<double> zA_plot(zA.begin(), zA.begin() + plot_N);
    std::vector<double> zP_plot(zP.begin(), zP.begin() + plot_N);
    std::vector<double> zF_plot(zF.begin(), zF.begin() + plot_N);


    plt::subplot(3, 1, 1);
    plt::plot(t_plot, zA_plot, "gray");
    plt::ylabel("Amplituda");
    plt::xlim(0.0, 0.30);
    plt::ylim(-1.3, 1.4);

    plt::subplot(3, 1, 2);
    plt::plot(t_plot, zP_plot, "gray");
    plt::ylabel("Amplituda");
    plt::xlim(0.0, 0.30);
    plt::ylim(-1.3, 1.4);

    plt::subplot(3, 1, 3);
    plt::plot(t_plot, zF_plot, "gray");
    plt::xlabel("Czas [s]");
    plt::ylabel("Amplituda");
    plt::xlim(0.0, 0.30);
    plt::ylim(-1.3, 1.4);

    plt::save("zadanie_1.png");
    plt::show();

    return 0;
}