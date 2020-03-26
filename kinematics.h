#ifndef __KINEMATICS_H
#define __KINEMATICS_H

#include "linear_algebra.h"
#include "util.h"
#include "transformations.h"

namespace rkd {
	vector<vector<float>> jacobian(vector<vector<float>> args);

	vector<float> z(int i, vector<vector<vector<float>>> Ts);

	vector<float> p(int i, vector<vector<vector<float>>> Ts);
}

#endif // __KINEMATICS_H
