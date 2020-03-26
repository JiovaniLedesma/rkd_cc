#include "linear_algebra.h"
#include "util.h"

using namespace rkd;

void rkd::sv_mult(float scalar, vector<float> vec, vector<float> &vec_out)
{
	for (int i = 0; i < vec.size(); ++i)
		vec[i] *= scalar;

	vec_out = vec;
}


vector<float> rkd::sv_mult(float scalar, vector<float> vec)
{
	for (int i = 0; i < vec.size(); ++i)
		vec[i] *= scalar;

	return vec;
}


void rkd::m_dot(vector<vector<float>> mat1, vector<vector<float>> mat2, vector<vector<float>> &mat_out)
{
	int nRow1 = mat1.size(), nCol1 = mat1[0].size(), nRow2 = mat2.size(), nCol2 = mat2[0].size();

	if ((nCol1 == nRow2))
	{
		vector<vector<float>> mat(nRow1, vector<float>(nCol2));

		for (int i = 0; i < nRow1; ++i)
		{
			for (int j = 0; j < nCol2; ++j)
			{
				mat_out[i][j] = 0;
				for (int k = 0; k < nCol1; ++k)
					mat[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
		mat_out = mat;
	}
	else
		print("El número de columnas de la primera matriz no es el mismo número de renglones de la segunda matriz");
}


vector<vector<float>> rkd::m_dot(vector<vector<float>> mat1, vector<vector<float>> mat2)
{
	int nRow1 = mat1.size(), nCol1 = mat1[0].size(), nRow2 = mat2.size(), nCol2 = mat2[0].size();

	vector<vector<float>> mat_out(nRow1, vector<float>(nCol2));

	if ((nCol1 == nRow2))
	{
		vector<vector<float> >mat(nRow1, vector<float>(nCol2));

		for (int i = 0; i < nRow1; ++i)
		{
			for (int j = 0; j < nCol2; ++j)
			{
				mat_out[i][j] = 0;
				for (int k = 0; k < nCol1; ++k)
					mat[i][j] += mat1[i][k] * mat2[k][j];
			}
		}

		mat_out = mat;
	}
	else
		print("El número de columnas de la primera matriz no es el mismo número de renglones de la segunda matriz");

	return mat_out;
}

void rkd::m_mult(vector<vector<vector<float>>> matrices, vector<vector<float>> &mat_out)
{
	int f = matrices[0].size();
	mat_out = eye(f);

	for (int i = 0; i < matrices.size(); ++i)
		mat_out = m_dot(mat_out, matrices[i]);
}

vector<vector<float>> rkd::m_mult(vector<vector<vector<float>>> matrices)
{
	int f = matrices[0].size();

	vector<vector<float>> mat = eye(f);

	for (int i = 0; i < matrices.size(); ++i)
		mat = m_dot(mat, matrices[i]);

	return mat;
}

void rkd::transpose(vector<vector<float>> mat, vector<vector<float>> &mat_out)
{
	int nRow = mat.size(), nCol = mat[0].size();

	vector<vector<float>> trans_mat(nCol, vector<float>(nRow));

	for (int i = 0; i < nRow; ++i)
		for (int j = 0; j < nCol; ++j)
			trans_mat[j][i] = mat[i][j];

	mat_out = trans_mat;
}

vector<vector<float>> rkd::transpose(vector<vector<float>> mat)
{
	int nRow = mat.size(), nCol = mat[0].size();

	vector<vector<float>> trans_mat(nCol, vector<float>(nRow));

	for (int i = 0; i < nRow; ++i)
		for (int j = 0; j < nCol; ++j)
			trans_mat[j][i] = mat[i][j];
	
	return trans_mat;
}

void rkd::zeros(int ndim, vector<vector<float>> &mat_out)
{
	vector<float> vec;

	for (int i = 0; i < ndim; ++i){
		for (int j = 0; j < ndim; ++j)
			vec.push_back(0);

		mat_out.push_back(vec);
		vec.erase(vec.begin(), vec.end());
	}
}

void rkd::zeros(int m, int n, vector<vector<float>> &mat_out)
{
	vector<float> vec;

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j)
			vec.push_back(0);
		
		mat_out.push_back(vec);
		vec.erase(vec.begin(), vec.end());
	}

}

vector<vector<float>> rkd::zeros(int ndim)
{
	vector<vector<float>> mat(ndim, vector<float>(ndim));

	for (int i = 0; i < ndim; ++i)
		for (int j = 0; j < ndim; ++j)
			mat[i][j] = 0;

	return mat;
}

vector<vector<float>> rkd::zeros(int m, int n)
{
	vector<vector<float>> mat(m, vector<float>(n));

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = 0;

	return mat;
}

void rkd::ones(int ndim, vector<vector<float>> &mat_out)
{
	vector<float> vec;

	for (int i = 0; i < ndim; ++i) {
		for (int j = 0; j < ndim; ++j)
			vec.push_back(1);
		
		mat_out.push_back(vec);
		vec.erase(vec.begin(), vec.end());
	}
}

void rkd::ones(int m, int n, vector<vector<float>> &mat_out)
{
	vector<float> vec;

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j)
			vec.push_back(1);

		mat_out.push_back(vec);
		vec.erase(vec.begin(), vec.end());
	}
}

vector<vector<float>> rkd::ones(int ndim)
{
	vector<vector<float>> mat(ndim, vector<float>(ndim));

	for (int i = 0; i < ndim; ++i)
		for (int j = 0; j < ndim; ++j)
			mat[i][j] = 1;
	
	return mat;
}

vector<vector<float>> rkd::ones(int m, int n)
{
	vector<vector<float>> mat(m, vector<float>(n));

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = 1;

	return mat;
}

void rkd::eye(int ndim, vector<vector<float>> &mat_out)
{
	mat_out = zeros(ndim);

	for (int i = 0; i < ndim; ++i)
		for (int j = 0; j < ndim; ++j)
			if (i == j)
				mat_out[i][j] = 1;
}

vector<vector<float>> rkd::eye(int ndim)
{
	vector<vector<float>> I = zeros(ndim);

	for (int i = 0; i < ndim; ++i)
		for (int j = 0; j < ndim; ++j)
			if (i == j)
				I[i][j] = 1;
	return I;
}

void rkd::cross(vector<float> a, vector<float> b, vector<float> &vec_out)
{
	float i = ((a[1] * b[2]) - (a[2] - b[1]));
	float j = -1 * ((a[0] * b[2]) - (a[2] - b[0]));
	float k = ((a[0] * b[1]) - (a[1] - b[0]));
	
	vec_out = { i, j, k };
}

vector<float> rkd::cross(vector<float> a, vector<float> b)
{
	float i = ((a[1] * b[2]) - (a[2] * b[1]));
	float j = -1 * ((a[0] * b[2]) - (a[2] * b[0]));
	float k = ((a[0] * b[1]) - (a[1] * b[0]));

	vector<float> vec_out = { i, j, k };

	return vec_out;
}

void rkd::roi_mat(vector<vector<float>> mat_in, vector<vector<float>> &mat_out, int x1, int x2, int y1, int y2)
{
	int subst1 = x2 - x1;
	int subst2 = y2 - y1;
	
	vector<float> vec;

	if ((subst1 > 0) && (subst2 > 0))
	{
		for (int i = y1; i < y2; ++i) {
			for (int j = x1; j < x2; ++j)
				vec.push_back(mat_in[i][j]);
			
			mat_out.push_back(vec);
			vec.erase(vec.begin(), vec.end());
		}
	}
	else
		print("Las variables de x1, x2, y1, y2; no deben de tener el mismo valor.\nEs decir, x1 < x2; y1 < y2");
}

void rkd::roi_mat2vec(vector<vector<float>> mat, vector<float> &vec_out, int x1, int x2, int y1, int y2)
{
	int subst1 = x2 - x1;
	int subst2 = y2 - y1;

	if ((subst1 == 1) || (subst2 == 1))
	{
		for (int i = y1; i < y2; ++i)
			for (int j = x1; j < x2; ++j)
				vec_out.push_back(mat[i][j]);
	}
}

vector<vector<float>> rkd::roi_mat(vector<vector<float>> mat, int x1, int x2, int y1, int y2)
{
	int subst1 = x2 - x1;
	int subst2 = y2 - y1;

	vector<vector<float>> mat_out;
	vector<float> vec;

	if ((subst1 > 0) && (subst2 > 0))
	{
		for (int i = y1; i < y2; ++i) {
			for (int j = x1; j < x2; ++j)
				vec.push_back(mat[i][j]);
			
			mat_out.push_back(vec);
			vec.erase(vec.begin(), vec.end());
		}
	}

	return mat_out;
}

vector<float> rkd::roi_mat2vec(vector<vector<float>> mat, int x1, int x2, int y1, int y2)
{
	int subst1 = x2 - x1;
	int subst2 = y2 - y1;

	vector<float> vec;

	if ((subst1 == 1) || (subst2 == 1))
	{
		for (int i = y1; i < y2; ++i)
			for (int j = x1; j < x2; ++j)
				vec.push_back(mat[i][j]);
	}

	return vec;
}

void rkd::reshape(vector<vector<float>> mat, vector<vector<float>> &mat_out, int m, int n)
{
	int nRow = mat.size(), nCol = mat[0].size();
	
	vector<vector<float>> Mat(m, vector<float>(n));

	if (nRow * nCol == m * n)
	{
		int k = 0, l = 0;
		for (int i = 0; i < nRow; ++i)
		{
			for (int j = 0; j < nCol; ++j)
			{
				Mat[k][l++] = mat[i][j];
				if (j >= n)
				{
					k++;
					l = 0;
				}
			}
		}
	}

	mat_out = Mat;
}

void rkd::reshape(vector<float> vec, vector<vector<float>> &mat_out, int nCols)
{
	int elements = vec.size();
	
	const auto nRows = elements / nCols;

	const auto begin = std::begin(vec);

	if (elements == nRows * nCols)
	{
		for (std::size_t row = 0; row < nRows; ++row)
		{
			mat_out.push_back( {  begin + row*nCols, begin + (row+1)*nCols } );
		}
	}
}

vector<vector<float>> rkd::reshape(vector<vector<float>> mat, int m, int n)
{
	int nRow = mat.size(), nCol = mat[0].size();

	vector<vector<float>> Mat(nRow, vector<float>(nCol));

	if (nRow * nCol == m * n)
	{
		int k = 0, l = 0;

		for (int i = 0; i < nRow; ++i)
		{
			for (int j = 0; j < nCol; ++j)
			{
				Mat[k][l++] = mat[i][j];
				if (j >= n)
				{
					k++;
					l = 0;
				}
			}
		}
	}

	return Mat;
}

vector<vector<float>> rkd::reshape(vector<float> vec, int nCols)
{
	int elements = vec.size();

	const auto nRows = elements / nCols;

	vector<vector<float>> mat;

	const auto begin = std::begin(vec);
	
	if (elements == nRows * nCols)
	{
		for (std::size_t row = 0; row < nRows; ++row)
		{
			mat.push_back( { begin + row*nCols, begin + (row+1)*nCols } );
		}
	}

	return mat;
}

void rkd::row_stack(vector<float> a, vector<float> b, vector<vector<float>> &mat_out)
{
	mat_out.push_back( a );
	mat_out.push_back( b );
}

void rkd::row_stack(vector<vector<float>> a, vector<vector<float>> b, vector<vector<float>> &mat_out)
{
	int nRow1 = a.size(), nCol1 = a[0].size(), nRow2 = b.size(), nCol2 = b[0].size();

	vector<float> vec;

	if ( nRow1 >= 1 && nRow2 >= 1)
	{
		if ( nCol1 == 1 && nCol2 == 1)
		{

			for (int i = 0; i < nRow1; ++i) {
				for (int j = 0; j < nCol1; ++j)
					vec.push_back(a[i][j]);
				
				mat_out.push_back(vec);
				vec.erase(vec.begin(), vec.end());
			}

			for (int i = 0; i < nRow2; ++i) {
				for (int j = 0; j < nCol2; ++j)
					vec.push_back(b[i][j]);

				mat_out.push_back(vec);
				vec.erase(vec.begin(), vec.end());
			}
		}
		else
		{
			for (int i = 0; i < nRow1; ++i)
				for (int j = 0; j < nCol1; ++j)
					vec.push_back(a[i][j]);

			mat_out.push_back( vec );

			vec.erase(vec.begin(), vec.end());

			for (int i = 0; i < nRow2; ++i)
				for (int j = 0; j < nCol2; ++j)
					vec.push_back(b[i][j]);
			
			mat_out.push_back( vec );
		}
	}
}

vector<vector<float>> rkd::row_stack(vector<float> a, vector<float> b)
{
	vector<vector<float>> mat;

	mat.push_back( a );
	mat.push_back( b );

	return mat;
}

vector<vector<float>> rkd::row_stack(vector<vector<float>> a, vector<vector<float>> b)
{
	int nRow1 = a.size(), nCol1 = a[0].size(), nRow2 = b.size(), nCol2 = b[0].size();
 
	vector<vector<float>> mat;
	vector<float> vec;

    if ( nRow1 >= 1 && nRow2 >= 1)
    {
        if ( nCol1 == 1 && nCol2 == 1)
        {
            for (int i = 0; i < nRow1; ++i) {
                for (int j = 0; j < nCol1; ++j)
                    vec.push_back(a[i][j]);

                mat.push_back(vec);
                vec.erase(vec.begin(), vec.end());
            }
 
            for (int i = 0; i < nRow2; ++i) {
                for (int j = 0; j < nCol2; ++j)
                    vec.push_back(b[i][j]);
 
                mat.push_back(vec);
                vec.erase(vec.begin(), vec.end());
            }
        }
        else
        {
            for (int i = 0; i < nRow1; ++i)
                for (int j = 0; j < nCol1; ++j)
                   vec.push_back(a[i][j]);

            mat.push_back( vec );

			vec.erase(vec.begin(), vec.end());

            for (int i = 0; i < nRow2; ++i)
                for (int j = 0; j < nCol2; ++j)
                    vec.push_back(b[i][j]);

            mat.push_back( vec );
		}
	}

	return mat;
}

void rkd::subtractV(vector<float> a, vector<float> b, vector<float> &vec_out)
{
	int na = a.size(), nb = b.size();

	if (na == nb)
	{
		for (int i = 0; i < na; ++i)
			vec_out.push_back(a[i] - b[i]);
	}
}

vector<float> rkd::subtractV(vector<float> a, vector<float> b)
{
	int na = a.size(), nb = b.size();
	vector<float> vec;

	if (na == nb)
	{
		for (int i = 0; i < na; ++i)
			vec.push_back(a[i] - b[i]);
	}

	return vec;
}

void rkd::flatten(vector<vector<float>> mat, vector<float> &vec_out)
{
	int nRow = mat.size(), nCol = mat[0].size();

	for (int i = 0; i < nRow; ++i)
		for (int j = 0; j < nCol; ++j)
			vec_out.push_back(mat[i][j]);
}

vector<float> rkd::flatten(vector<vector<float>> mat)
{
	int nRow = mat.size(), nCol = mat[0].size();
	
	vector<float> vec;

	for (int i = 0; i < nRow; ++i)
		for (int j = 0; j < nCol; ++j)
			vec.push_back(mat[i][j]);

	return vec;
}

void rkd::mat2vec(vector<vector<float>> mat, vector<float> &vec_out)
{
	int nRow = mat.size(), nCol = mat[0].size();

	for (int i = 0; i < nRow; ++i)
		for (int j = 0; j < nCol; ++j)
			vec_out.push_back(mat[i][j]);
}

vector<float> rkd::mat2vec(vector<vector<float>> mat)
{
	int nRow = mat.size(), nCol = mat[0].size();

	vector<float> vec;

	for (int i = 0; i < nRow; ++i)
		for (int j = 0; j < nCol; ++j)
			vec.push_back(mat[i][j]);

	return vec;
}

void rkd::add_vec2mat(vector<float> vec, vector<vector<float>> mat, vector<vector<float>> &mat_out, bool row, bool col, int index)
{
	int nRow = mat.size(), nCol = mat[0].size(), nElem = vec.size();

	if (row)
	{
		if (nRow == nElem)
		{
			for (int i = 0; i < nRow; ++i)
				mat[i][index] = vec[i];
		}
	}
	else if (col)
	{
		if (nCol == nElem)
		{
			for (int i = 0; i < nCol; ++i)
				mat[index][i] = vec[i];
		}
	}

	mat_out = mat;
}

vector<vector<float>> rkd::add_vec2mat(vector<float> vec, vector<vector<float>> mat, bool row, bool col, int index)
{
	int nRow = mat.size(), nCol = mat[0].size(), nElem = vec.size();

	if (row)
	{
		if (nRow == nElem)
		{
			for (int i = 0; i < nRow; ++i)
				mat[i][index] = vec[i];
		}
	}
	else if (col)
	{
		if (nCol == nElem)
		{
			for (int i = 0; i < nCol; ++i)
				mat[index][i] = vec[i];
		}
	}

	return mat;
}
