// table 3(6)
// table 4(5)
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

int main() {
    const double fs_u = 8000.0;
    const double Tc_u = 4.5;
    const int N_u = static_cast<int>(fs_u * Tc_u);

    vector<double> t_u(N_u), u(N_u);

    for (int i = 0; i < N_u; ++i) {
        t_u[i] = i / fs_u;
        double t = t_u[i];

        // table 3(6)
        if (t >= 0 && t < 1.8)           
            u[i] = -(t / 2.0) * sin(20.0 * pow(t, 3) - 18.0 * pow(t, 2));
        else if (t >= 1.8 && t < 3.0)    
            u[i] = cos(5.0 * M_PI * t) * sin(12.0 * M_PI * pow(t, 2));
        else if (t >= 3.0 && t < 4.5)    
            u[i] = ((t - 3.0) / 3.0) * sin((12.0 - t) * M_PI * pow(t, 2));
    }

    // table 4 (bk(t)) 
    const double fs_b = 22050.0;
    const double Tc_b = 1.0;
    const int N_b = static_cast<int>(fs_b * Tc_b);

    vector<int> H = { 2, 20, 40 };
    vector<double> t_b(N_b), b1(N_b), b2(N_b), b3(N_b);
    vector<double*> signals = { &b1[0], &b2[0], &b3[0] };

    for (int i = 0; i < N_b; ++i) {
        t_b[i] = i / fs_b;
        double t = t_b[i];

        // table 4(5)
        for (int k = 0; k < 3; ++k) {
            double sum = 0;
            for (int h = 1; h <= H[k]; ++h) {
                sum += (pow(-1, h) / (3.0 * pow(h, 2))) * cos(2.0 * M_PI * h * t + sin(6.0 * M_PI * t));
            }
            signals[k][i] = sum;
        }
    }

    plt::figure();

    plt::subplot(2, 2, 1); plt::plot(t_u, u); plt::title("Tab. 3, Var. 6: u(t)"); plt::grid(true);
    plt::subplot(2, 2, 2); plt::plot(t_b, b1); plt::title("Tab. 4, Var. 5: b1(t)"); plt::grid(true);
    plt::subplot(2, 2, 3); plt::plot(t_b, b2); plt::title("Tab. 4, Var. 5: b2(t)"); plt::grid(true);
    plt::subplot(2, 2, 4); plt::plot(t_b, b3); plt::title("Tab. 4, Var. 5: b3(t)"); plt::grid(true);

    plt::save("zadanie_2.png");
    plt::show();
    return 0;
}