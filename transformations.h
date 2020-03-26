#ifndef __TRANSFORMATIONS_H
#define __TRANSFORMATIONS_H

#include <vector>
#include <cmath>

using std::vector;

namespace rkd{
vector<vector<float>> rotx(float theta, bool deg=false);
vector<vector<float>> roty(float theta, bool deg=false);
vector<vector<float>> rotz(float theta, bool deg=false);

vector<float> rot2eul(vector<vector<float>> mat, const char *axis, bool = false, bool sol = false);
vector<float> rot2RPY(vector<vector<float>> mat, bool deg = false, bool sol = false);
vector<vector<float>> rot2axa(vector<vector<float>> mat, bool deg = false);
vector<vector<float>> RPY2Rot(float phi, float theta, float psi, bool deg = false);

vector<vector<float>> htmDH(float a, float alpha, float d, float t, bool deg = false);
}

#endif // __TRANSFORMATIONS_H
