#ifndef __UTIL_H
#define __UTIL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include "TextTable.h"

#define pi M_PI

#define alpha_symbol     "\u03B1 "
#define beta_symbol      "\u03B2 "
#define gamma_symbol     "\u03B3 "
#define delta_symbol     "\u03B4 "
#define epsilon_symbol   "\u03B5 "
#define zeta_symbol      "\u03B6 "
#define eta_symbol       "\u03B7 "
#define theta_symbol     "\u03B8 "
#define iota_symbol      "\u03B9 "
#define kappa_symbol     "\u03BA "
#define lamda_symbol     "\u03BB "
#define mu_symbol        "\u03BC "
#define nu_symbol        "\u03BD "
#define xi_symbol        "\u03BE "
#define omicron_symbol   "\u03BF "
#define pi_symbol        "\u03C1 "
#define rho_symbol       "\u03C2 "
#define sigma_symbol     "\u03C3 "
#define tau_symbol       "\u03C4 "
#define upsilon_symbol   "\u03C5 "
#define phi_symbol       "\u03C6 "
#define chi_symbol       "\u03C7 "
#define psi_symbol       "\u03C8 "
#define omega_symbol     "\u03C9 "

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::to_string;

namespace rkd{

void deg2rad(float &x);
void rad2deg(float &x);
void print(vector<float> vec);
void print(vector<vector<float>> mat, vector<string> titles);
void print(vector<vector<float>> mat);
void print(const char *cadena);
void print(int value);
void print(float value);
void print(double value);
}

#endif // __UTIL_H
