#include "../vlfeat.hpp"
#include <iostream>

using namespace arma;
using namespace std;

//void vl_kmeans(const arma::fmat &X, arma::fmat &Y, // for assignments and energy
//        int numCenters,
//        const std::string &algorithm_ = "lloyd",
//        const std::string &distance_ = "l2",
//        const std::string &initialization_ = "plusplus",
//        double minEnergyVariation = -1,
//        int numRepetitions = 1,
//        int numTrees = 3,
//        int maxNumComparisons = 100,
//        int maxNumIterations = 100,
//        int verbose = 0)

int main() {
	fmat X, Y;
	
	X = arma::randn<fmat>(5,100);
	Y = arma::fmat();
	
	vl_kmeans(X, Y, 10);

	cout << "X shape: " << X.n_rows << " "<< X.n_cols << endl;
	cout << "Y shape: " << Y.n_rows << " "<< Y.n_cols << endl;

	return 0;
}
