#include "transformations.h"
#include "util.h"
#include "linear_algebra.h"

using namespace rkd;

vector<vector<float>> rkd::rotx(float theta, bool deg)
{
    if (deg)
        rkd::deg2rad(theta);
    
    float ct = cos(theta);
    float st = sin(theta);

    vector<vector<float>> mat = { {1,  0,   0},
                                  {0, ct, -st},
                                  {0, st,  ct}  };
    
    return mat;
}

vector<vector<float>> rkd::roty(float theta, bool deg)
{
    if (deg)
        rkd::deg2rad(theta);
    
    float ct = cos(theta);
    float st = sin(theta);

    vector<vector<float>> mat = { {ct,  0, st},
                                  { 0,  1,  0},
                                  {-st, 0, ct}  };
    
    return mat;
}

vector<vector<float>> rkd::rotz(float theta, bool deg)
{
    if (deg)
        rkd::deg2rad(theta);
    
    float ct = cos(theta);
    float st = sin(theta);

    vector<vector<float>> mat = { {ct, -st, 0},
                                  {st,  ct, 0},
                                  { 0,   0, 1}  };
    
    return mat;
}

vector<float> rkd::rot2eul(vector<vector<float>> mat, const char *axis, bool deg, bool sol)
{
	float r11 = mat[0][0];
	float r12 = mat[0][1];
	float r13 = mat[0][2];

	float r21 = mat[1][0];
	float r22 = mat[1][1];
	float r23 = mat[1][2];

	float r31 = mat[2][0];
	float r32 = mat[2][1];
	float r33 = mat[2][2];

	float theta, phi, psi;
	
	vector<float> result;

	if ((axis == "ZXZ") || (axis == "zxz"))
	{
		if ((r33 != 1) || (r33 != -1))
		{
			theta = atan2(sqrt(1-(pow(r33, 2))), r33);
			phi   = atan2(r13, -r23);
			psi   = atan2(r31, r32);
		}

		if (sol) // This condition is the user want the second result
		{
			theta = atan2(-sqrt(1-pow(r33, 2)), r33);
			phi   = atan2(-r13, r23);
			psi   = atan2(-r31, -r32);
		}

		if (r33 == 1)
		{
			theta = 0;
			phi   = 0;
			psi   = atan2(r21, r11);
		}

		if (r33 == -1)
		{
			theta = pi;
			phi   = 0;
			psi   = atan2(-r21, -r11);
		}
	}

	if ((axis == "ZYZ") || (axis == "zyz"))
	{
		if ((r33 != 1) || (r33 != -1))
		{
			theta = atan2(sqrt(1-pow(r33, 2)), r33);
			phi   = atan2(r23, r13);
			psi   = atan2(r32, -r31);
		}

		if (sol)
		{
			theta = atan2(-sqrt(1-pow(r33, 2)), r33);
			phi   = atan2(-r23, -r13);
			psi   = atan2(-r32, r31);
		}

		if (r33 == 1)
		{
			theta = 0;
			phi   = 0;
			psi   = atan2(r21, r11);
		}

		if (r33 == -1)
		{
			theta = pi;
			phi   = 0;
			psi   = atan2(-r21, -r11);
		}
	}

	if ((axis == "XZX") || (axis == "xzx"))
	{
		if ((r11 != 1) || (r11 != -1))
		{
			theta = atan2(sqrt(1-pow(r11, 2)), r11);
			phi   = atan2(r31, r21);
			psi   = atan2(r13, -r12);
		}

		if (sol)
		{
			theta = atan2(-sqrt(1-pow(r11, 2)), r11);
			phi   = atan2(-r31, -r21);
			psi   = atan2(-r13, r12);
		}

		if (r11 == 1)
		{
			theta = 0;
			phi   = 0;
			psi   = atan2(r32, r22);
		}

		if (r11 == -1)
		{
			theta = pi;
			phi   = 0;
			psi   = atan2(-r32, r22);
		}
	}

	if ((axis == "XYX") || (axis == "xyx"))
	{
		if ((r11 != 1) || (r11 != -1))
		{
			theta = atan2(sqrt(1-pow(r11, 2)), r11);
			phi   = atan2(r21, -r31);
			psi   = atan2(r12, r13);
		}

		if (sol)
		{
			theta = atan2(-sqrt(1-pow(r11, 2)), r11);
			phi   = atan2(-r21, r31);
			psi   = atan2(-r12, r13);
		}

		if (r11 == 1)
		{
			theta = 0;
			phi   = 0;
			psi   = atan2(r32, r22);
		}

		if (r11 == -1)
		{
			theta = pi;
			phi   = 0;
			psi   = atan2(-r32, -r22);
		}
	}

	if ((axis == "YZY") || (axis == "yzy"))
	{
		if ((r22 != 1) || (r22 != -1))
		{
			theta = atan2(sqrt(1-pow(r22, 2)), r22);
			phi   = atan2(r32, -r12);
			psi   = atan2(r23, r21);
		}

		if (sol)
		{
			theta = atan2(-sqrt(1-pow(r22, 2)), r22);
			phi   = atan2(-r32, r12);
			psi   = atan2(-r23, -r21);
		}

		if (r22 == 1)
		{
			theta = 0;
			phi   = 0;
			psi   = atan2(r13, r33);
		}

		if (r22 == -1)
		{
			theta = pi;
			phi   = 0;
			psi   = atan2(-r13, -r33);
		}
	}

	if ((axis == "YXY") || (axis == "yxy"))
	{
		if ((r22 != 1) || (r22 != -1))
		{
			theta = atan2(sqrt(1-pow(r22, 2)), r22);
			phi   = atan2(r12, r32);
			psi   = atan2(r21, -r23);
		}

		if (sol)
		{
			theta = atan2(-sqrt(1-pow(r22, 2)), r22);
			phi   = atan2(-r12, r32);
			psi   = atan2(-r21, -r23);
		}

		if (r22 == 1)
		{
			theta = 0;
			phi   = 0;
			psi   = atan2(r13, r33);
		}

		if (r22 == -1)
		{
			theta = pi;
			phi   = 0;
			psi   = atan2(-r13, -r33);
		}
	}

	if (deg)
	{
		rad2deg(phi);
		rad2deg(theta);
		rad2deg(psi);
	}
	
	result.push_back(phi);
	result.push_back(theta);
	result.push_back(psi);

	return result;
}

vector<float> rkd::rot2RPY(vector<vector<float>> mat, bool deg, bool sol)
{
	float r11 = mat[0][0];
    float r12 = mat[0][1];
    float r13 = mat[0][2];
  
	float r21 = mat[1][0];
    float r22 = mat[1][1];
    float r23 = mat[1][2];
  
    float r31 = mat[2][0];
    float r32 = mat[2][1];
    float r33 = mat[2][2];
 
    float theta, phi, psi;
  
    vector<float> result;

	if ((r31 != 1) || (r31 != -1))
	{
		theta = atan2(r31, sqrt(1-pow(r31, 2)));
		phi   = atan2(r21, r11);
		psi   = atan2(r32, r33);
	}

	if (sol)
	{
		theta = atan2(r31, -sqrt(1-pow(r31, 2)));
		phi   = atan2(-r21, -r11);
		psi   = atan2(-r31, -r33);
	}

	if (r31 == 1)
	{
		theta = pi/2;
		phi   = 0;
		psi   = atan2(-r21, r22);
	}

	if (r31 == -1)
	{
		theta = (3*pi)/2;
		phi   = 0;
		psi   = atan2(r21, -r22);
	}

	if (deg)
	{
		rad2deg(phi);
		rad2deg(theta);
		rad2deg(psi);
	}

	result.push_back(phi);
	result.push_back(theta);
	result.push_back(psi);

	return result;
}

vector<vector<float>> rkd::rot2axa(vector<vector<float>> mat, bool deg)
{
	float r11 = mat[0][0];
    float r12 = mat[0][1];
    float r13 = mat[0][2];
    float r21 = mat[1][0];
    float r22 = mat[1][1];
    float r23 = mat[1][2]; 
    float r31 = mat[2][0];
    float r32 = mat[2][1];
    float r33 = mat[2][2];

	std::cout << "r11: " << r11 << "\tr12: " << r12 << "\tr13: " << r13 << std::endl;
	std::cout << "r21: " << r21 << "\tr12: " << r22 << "\tr33: " << r23 << std::endl;
	std::cout << "r31: " << r31 << "\tr12: " << r32 << "\tr13: " << r33 << std::endl;
  
    float theta;
   
    vector<float> angle;
	vector<vector<float>> result;

	theta = acos((r11 + r22 + r33 - 1) / 2);
	
	float sc = 1/(2*sin(theta));
	vector<float> vec = { r32-r23, r13-r31, r21-r12 };
	
	vector<float> k = sv_mult(sc, vec);

	if (deg)
		rad2deg(theta);
	
	angle.push_back(theta);

	result.push_back(angle);
	result.push_back(k);

	return result;
}

vector<vector<float>> rkd::RPY2Rot(float phi, float theta, float psi, bool deg)
{
	vector<vector<float>> Rz = rotz(phi, deg);
	vector<vector<float>> Ry = roty(theta, deg);
	vector<vector<float>> Rx = rotx(psi, deg);
	
	vector<vector<vector<float>>> matrices = { Rz, Ry, Rx };
	vector<vector<float>> mat;
	m_mult(matrices, mat);

	return mat;
}

vector<vector<float>> rkd::htmDH(float a, float alpha, float d, float t, bool deg)
{
	if (deg)
	{
		deg2rad(alpha);
		deg2rad(t);
	}

	float cal = cos(alpha);
	float sal = sin(alpha);
	float ct  = cos(t);
	float st  = sin(t);

	vector<vector<float>> H = { { ct, -st*cal,  st*sal, a*ct },
   								{ st,  ct*cal, -ct*sal, a*st },
								{ 0,      sal,     cal,    d },
								{ 0,        0,       0,    1 }	};

	return H;
}
