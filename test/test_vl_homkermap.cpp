#include "../vlfeat.hpp"
#include <iostream>

using namespace arma;
using namespace std;

//void vl_homkermap(const arma::fmat &XX, arma::fmat &VV, int n,
//        const std::string & kernel = "kchi2",
//        const std::string & window = "rectangular",
//        double gamma = 1.0,
//        double period = -1);

int main() {
	fmat X, V;
	
	X = arma::randn<fmat>(3,10);
	vl_homkermap(X, V, 1);

	cout << "X shape: " << X.n_rows << " "<< X.n_cols << endl;
	cout << "V shape: " << V.n_rows << " "<< V.n_cols << endl;
	
	cout << X << endl;
	cout << V << endl;

	return 0;
}
