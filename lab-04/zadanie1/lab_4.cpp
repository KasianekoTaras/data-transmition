#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    double fm = 2.0;
    double fn = 50.0;
    double fs = 1000.0;
    double duration = 1.0;

    double dt = 1.0 / fs;
    int num_samples = duration * fs;

    std::vector<double> t(num_samples);
    std::vector<double> m(num_samples);

    for (int i = 0; i < num_samples; ++i) {
        t[i] = i * dt;
        m[i] = std::sin(2 * M_PI * fm * t[i]);
    }

    int plot_idx = 1;

    std::vector<double> kA_vals = { 0.5, 5.0, 25.0 };
    for (double kA : kA_vals) {
        std::vector<double> zA(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            zA[i] = (kA * m[i] + 1.0) * std::cos(2 * M_PI * fn * t[i]);
        }
        plt::subplot(3, 3, plot_idx++);
        plt::plot(t, zA, "gray");
        plt::plot(t, m, "blue");
        plt::title("AM: kA = " + std::to_string(kA));
    }

    std::vector<double> kP_vals = { 0.5, 2.0, 8.0 };
    for (double kP : kP_vals) {
        std::vector<double> zP(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            zP[i] = std::cos(2 * M_PI * fn * t[i] + kP * m[i]);
        }
        plt::subplot(3, 3, plot_idx++);
        plt::plot(t, zP, "gray");
        plt::plot(t, m, "blue");
        plt::title("PM: kP = " + std::to_string(kP));
    }

    std::vector<double> kF_vals = { 0.5, 2.0, 8.0 };
    for (double kF : kF_vals) {
        std::vector<double> zF(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            zF[i] = std::cos(2 * M_PI * fn * t[i] + (kF / fm) * m[i]);
        }
        plt::subplot(3, 3, plot_idx++);
        plt::plot(t, zF, "gray");
        plt::plot(t, m, "blue");
        plt::title("FM: kF = " + std::to_string(kF));
    }

    plt::tight_layout();
    plt::save("zadanie_1.png");
    plt::show();

    return 0;
}