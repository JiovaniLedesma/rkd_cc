#include <chrono> // To measure time
#include <vector>
#include <iostream>
#include "transformations.h"
#include "linear_algebra.h"
#include "kinematics.h"
#include "util.h"

using std::vector;
using namespace rkd;
using namespace std::chrono;

int main()
{
    auto start = high_resolution_clock::now(); 

	vector<float> DH1 = {0, pi/2, 100, 0.785398163, 1};
	vector<float> DH2 = {0, pi/2,   0, 1.047197551, 1};
	vector<float> DH3 = {0, 0, 150, 0, 0};
	vector<float> DH4 = {40, 50, 120, 0, 0};
	vector<float> DH5 = {0, 3, 45, 0, 1};
	vector<float> DH6 = {30, 25, 60, 3, 0};

	vector<vector<float>> result = jacobian( { DH1, DH2, DH3, DH4, DH5, DH6 } );

	print(">>> Jacobian Matrix <<<");
	print(result);

	vector<vector<float>> DH = htmDH(0, pi/2, 0.785398163, true);
	vector<string> titles = {"a", alpha_symbol, "d", theta_symbol};
	print("\n\n>>> DH MATRIX <<<");
	print(DH, titles);

    auto stop = high_resolution_clock::now();
	
	auto duration = duration_cast<microseconds>(stop - start);

    cout << "\n\n" << duration.count() << " ms" << endl;
}
