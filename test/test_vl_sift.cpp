#include "../vlfeat.hpp"
#include <iostream>

using namespace arma;
using namespace std;

//void vl_sift(const arma::fmat &I, arma::mat &f, arma::fmat &d,
//        int octaves = -1,
//        int levels = -1,
//        int first_octave = -1,
//        double peak_thresh = -1,
//        double edge_thresh = -1,
//        double norm_thresh = -1,
//        double magnif = -1,
//        double window_size = -1,
//        bool orientations = false,
//        bool float_descriptors = false,
//        int verbose = 0);
        
int main() {
	fmat I;
	mat f;
	fmat d;

	I.load("lena.txt", raw_ascii);
	cout << "image shape: " << I.n_rows << " "<< I.n_cols << endl;

//	vl_sift(I, f, d);
	vl_sift(I, f, d, -1, -1, -1, -1, -1, -1, -1, -1, false, true, 1);

	cout << "frames shape: " << f.n_rows << " "<< f.n_cols << endl;
	cout << "descr shape: " << d.n_rows << " "<< d.n_cols << endl;

	f.save("sift_f_cpp.txt", raw_ascii);
	d.save("sift_d_cpp.txt", raw_ascii);
	return 0;
}
