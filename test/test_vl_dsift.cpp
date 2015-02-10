#include "vlfeat.hpp"
#include <iostream>

using namespace arma;
using namespace std;

//void vl_dsift(const arma::fmat &I, arma::mat &f, arma::fmat &d,
//        const arma::mat &bounds_ = arma::mat(),
//        const arma::mat &step_ = arma::mat(),
//        const arma::mat &size_ = arma::mat(),
//        const arma::mat &geometry_ = arma::mat(),
//        bool fast = true,
//        bool norm = true,
//        double windowSize = -1,
//        bool floatDescriptors = false,
//        int verbose = 0);

int main() {
	fmat I;
	mat f;
	fmat d;

	I.load("lena.txt", raw_ascii);
	cout << "image shape: " << I.n_rows << " "<< I.n_cols << endl;

	double step[1] = {100};
	vl_dsift(I, f, d, mat(), mat(step,1,1), mat(), mat(), true, true, -1, true, 1);

	cout << "frames shape: " << f.n_rows << " "<< f.n_cols << endl;
	cout << "descr shape: " << d.n_rows << " "<< d.n_cols << endl;

	f.save("dsift_f_cpp.txt", raw_ascii);
	d.save("dsift_d_cpp.txt", raw_ascii);
	return 0;
}
