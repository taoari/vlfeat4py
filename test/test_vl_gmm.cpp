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

//void vl_gmm (const arma::fmat &X, int numClusters,
//	arma::fmat &means, arma::fmat &covariances, arma::fmat &priors, arma::fmat &posteriors,
//	const std::string &initialization_ = "rand", 
//	const arma::fmat & initMeans_ = arma::fmat(), 
//	const arma::fmat & initPriors_ = arma::fmat(),
//	const arma::fmat & initCovariances_ = arma::fmat(),
//	int maxNumIterations = 100,
//	int numRepetitions = 1,
//	const arma::mat &covarianceBound_ = arma::fmat(),
//	int verbose = 0);

int main() {
	fmat I;
	mat f;
	fmat d;

	I.load("lena.txt", raw_ascii);
	cout << "image shape: " << I.n_rows << " "<< I.n_cols << endl;

	double step[1] = {5};
	vl_dsift(I, f, d, mat(), mat(step,1,1), mat(), mat(), true, true, -1, true, 1);

	cout << "frames shape: " << f.n_rows << " "<< f.n_cols << endl;
	cout << "descr shape: " << d.n_rows << " "<< d.n_cols << endl;
	
	fmat means, covariances, priors, posteriors;
	
	vl_gmm(d, 30, means, covariances, priors, posteriors);
	
	fmat enc;
	vl_fisher(d, means, covariances, priors, enc);

	return 0;
}
