// table 1(8) x(t) = (1-t)*sin(2πft+phi)*cos(4πt)
// table 2(6) y(t) = 3t^0.3 * sin(20πt) - sin(f/7 * πt)
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

int main() {
    const double fs = 8000.0;
    const double Tc = 1.0;
    const double f = 5.0;
    const double phi = 0.0;
    const int N = static_cast<int>(fs * Tc);

    vector<double> t(N), x(N), y(N), z(N), v(N);
    double max_y = -99999;


    for (int i = 0; i < N; ++i) {
        t[i] = i / fs;
        double time = t[i];

        // table 1(8) x(t) = (1-t)*sin(2πft+phi)*cos(4πt)
        x[i] = (1.0 - time) * sin(2.0 * M_PI * f * time + phi) * cos(4.0 * M_PI * time);

        // table 2(6) y(t) = 3t^0.3 * sin(20πt) - sin(f/7 * πt)
        y[i] = 3.0 * pow(time, 0.3) * sin(20.0 * M_PI * time) - sin((f / 7.0) * M_PI * time);

        if (y[i] > max_y) max_y = y[i];
    }

    for (int i = 0; i < N; ++i) {
        double time = t[i];
        z[i] = -pow(time, 2.0) * sqrt(abs(max_y + y[i] - x[i] / 5.0));
        v[i] = -x[i] * (abs(y[i] * time) * exp(-x[i]));
    }

    plt::figure();
    plt::subplot(2, 2, 1); plt::plot(t, x); plt::title("Tab. 1, Var. 8: x(t)"); plt::grid(true);
    plt::subplot(2, 2, 2); plt::plot(t, y); plt::title("Tab. 2, Var. 6: y(t)"); plt::grid(true);
    plt::subplot(2, 2, 3); plt::plot(t, z); plt::title("Tab. 2, Var. 6: z(t)"); plt::grid(true);
    plt::subplot(2, 2, 4); plt::plot(t, v); plt::title("Tab. 2, Var. 6: v(t)"); plt::grid(true);

    plt::save("zadanie_1.png");
    plt::show();
    return 0;
}