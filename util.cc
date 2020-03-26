#include "util.h"

void rkd::deg2rad(float &x)
{
    x = (x*pi)/180;
}

void rkd::rad2deg(float &x)
{
    x = (x * 180) / pi;
}

void rkd::print(vector<float> vec)
{
    for (int i = 0; i < vec.size(); ++i)
        cout << vec[i] << endl;
}

void rkd::print(vector<vector<float>> mat, vector<string> titles)
{
    int nRow = mat.size(), nCol = mat[0].size();
	
	TextTable t( '-', '|', '+' );

	for (int i = 0; i < titles.size(); ++i)
	{
		t.add(titles[i]);
	}
	t.endOfRow();

    for (int i = 0; i < nRow; ++i)
    {
        for (int j = 0; j < nCol; ++j)
			t.add(to_string(mat[i][j]));
		t.endOfRow();
	}
	t.setAlignment( 2, TextTable::Alignment::LEFT );
    cout << t;
}

void rkd::print(vector<vector<float>> mat)
{
	int nRow = mat.size(), nCol = mat[0].size();
	
	TextTable t( '-', '|', '+' );

	for (int i = 0; i < nRow; ++i)
	{
		for (int j = 0; j < nCol; ++j)
			t.add(to_string(mat[i][j]));
		t.endOfRow();
	}
	t.setAlignment( 2, TextTable::Alignment::LEFT );
	cout << t;
}

void rkd::print(const char *cadena)
{
    cout << cadena << endl;
}

void rkd::print(int value)
{
    cout << value << endl;
}

void rkd::print(float value)
{
    cout << value << endl;
}

void rkd::print(double value)
{
    cout << value << endl;
}
