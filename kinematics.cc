#include "kinematics.h"

vector<vector<float>> rkd::jacobian(vector<vector<float>> args)
{
	vector<vector<vector<float>>> Ts;
	
	int dof = args.size();

	vector<float> k;
	vector<vector<float>> DH;

	int n = dof;
	vector<vector<float>> M_ = zeros(6, n);

	vector<float> jp, jo;
	vector<vector<float>> Jp;
	const char* type;

	for (int i = 0; i < dof; ++i)
	{
		k = args[i];
		Ts.push_back(htmDH(k[0], k[1], k[2], k[3]));
		k.erase(k.begin(), k.end());
	}

	for (int i = 0; i < dof; ++i)
	{
		k = args[i];
		if (k.size() > 4)
		{
			if (k[4] == 1)
				type = "r";
			else
				type = "p";
		}
		else
			type = "r";

		if ( type == "r" )
		{
			jp = cross(z(i, Ts), subtractV(p(n, Ts), p(i, Ts)));
			jo = z(i, Ts);
		}
		else
		{
			jp = z(i, Ts);
			jo = mat2vec(zeros(3, 1));
		}
		
		Jp = row_stack(mat2vec(reshape(jp, 1)), mat2vec(reshape(jo, 1)));
		jp.erase(jp.begin(), jp.end());
		jp = flatten(Jp);
		M_ = add_vec2mat(jp, M_, true, false, i);

		k.erase(k.begin(), k.end());
		jp.erase(jp.begin(), jp.end());
		jo.erase(jo.begin(), jo.end());
		Jp.erase(Jp.begin(), Jp.end());
	}

	return M_;
}

vector<float> rkd::z(int i, vector<vector<vector<float>>> Ts)
{
	if ( i == 0 )
		return { 0, 0, 1 };
	
	vector<vector<float>> MTH = eye(4);
	vector<vector<vector<float>>> matrices;

	for (int k = 0; k < i; ++k) {
		matrices = { MTH, Ts[k] };

		MTH = m_mult(matrices);

		matrices.erase(matrices.begin(), matrices.end());
	}
	
	vector<float> MTH_ = roi_mat2vec(MTH, 2, 3, 0, 3); // MTH[:3, 2]

	return MTH_;
}

vector<float> rkd::p(int i, vector<vector<vector<float>>> Ts)
{
	if (i == 0)
		return { 0, 0, 0 };
	
	vector<vector<float>> MTH = eye(4);
	vector<vector<vector<float>>> matrices;

	for (int k = 0; k < i; ++k) {
		matrices = { MTH, Ts[k] };

		MTH = m_mult(matrices);

		matrices.erase(matrices.begin(), matrices.end());
	}
	
	vector<float> MTH_ = roi_mat2vec(MTH, 3, 4, 0, 3); // MTH[:3, 3]

	return MTH_;
}
