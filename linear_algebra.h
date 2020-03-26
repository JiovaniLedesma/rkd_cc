#ifndef __LINERAR_ALGEBRA_H
#define __LINERAR_ALGEBRA_H

#include <vector>
#include <cmath>
#include <numeric>
#include "util.h"

using std::vector;


namespace rkd{
	void zeros(int ndim, vector<vector<float>> &mat_out);
	void zeros(int m, int n, vector<vector<float>> &mat_out);
	vector<vector<float>> zeros(int ndim);
	vector<vector<float>> zeros(int m, int n);
	
	void ones(int ndim, vector<vector<float>> &mat_out);
	void ones(int m, int n, vector<vector<float>> &mat_out);
	vector<vector<float>> ones(int ndim);
	vector<vector<float>> ones(int m, int n);
	
	void eye(int ndim, vector<vector<float>> &mat_out);
	vector<vector<float>> eye(int ndim);

	void cross(vector<float> a, vector<float> b, vector<float> &vec_out);
	vector<float> cross(vector<float> a, vector<float> b);

	void roi_mat(vector<vector<float>> mat_in, vector<vector<float>> &mat_out, int x1, int x2, int y1, int y2);
	void roi_mat2vec(vector<vector<float>> mat, vector<float> &vec_out, int x1, int x2, int y1, int y2);
	vector<vector<float>> roi_mat(vector<vector<float>> mat, int x1, int x2, int y1, int y2);
	vector<float> roi_mat2vec(vector<vector<float>> mat, int x1, int x2, int y1, int y2);

	void reshape(vector<vector<float>> mat, vector<vector<float>> &mat_out, int m, int n);
	void reshape(vector<float> vec, vector<vector<float>> &mat_out, int nCols);
	vector<vector<float>> reshape(vector<vector<float>> mat, int m, int n);
	vector<vector<float>> reshape(vector<float> vec, int nCols);

	void row_stack(vector<float> a, vector<float> b, vector<vector<float>> &mat_out);
	void row_stack(vector<vector<float>> a, vector<vector<float>> b, vector<vector<float>> &mat_out);
	vector<vector<float>> row_stack(vector<float> a, vector<float> b);
	vector<vector<float>> row_stack(vector<vector<float>> a, vector<vector<float>> b);

	void subtractV(vector<float> a, vector<float> b, vector<float> &vec_out);
	vector<float> subtractV(vector<float> a, vector<float> b);

	void flatten(vector<vector<float>> mat, vector<float> &vec_out);
	vector<float> flatten(vector<vector<float>> mat);

	void mat2vec(vector<vector<float>> mat, vector<float> &vec_out);
	vector<float> mat2vec(vector<vector<float>> mat);

	void add_vec2mat(vector<float> vec, vector<vector<float>> mat, vector<vector<float>> &mat_out, bool row, bool col, int index);
	vector<vector<float>> add_vec2mat(vector<float> vec, vector<vector<float>> mat, bool row, bool col, int index);

	void sv_mult(float scalar, vector<float> vec, vector<float> &vec_out);
	vector<float> sv_mult(float scalar, vector<float> vec);

	void m_dot(vector<vector<float>> mat1, vector<vector<float>> mat2, vector<vector<float>> &mat_out);
	vector<vector<float>> m_dot(vector<vector<float>> mat1, vector<vector<float>> mat2);

	void m_mult(vector<vector<vector<float>>> matrices, vector<vector<float>> &mat_out);
	vector<vector<float>> m_mult(vector<vector<vector<float>>> matrices);

	void transpose(vector<vector<float>> mat, vector<vector<float>> &mat_out);
	vector<vector<float>> transpose(vector<vector<float>> mat);
}

#endif // __LINERAR_ALGEBRA_H
