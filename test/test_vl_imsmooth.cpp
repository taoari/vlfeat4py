#include "../vlfeat.hpp"
#include <iostream>

using namespace arma;
using namespace std;

//void vl_imsmooth(const arma::mat &I, arma::mat &Is, 
//        double sigma;
//        const std::string &padding_ = "continuity",
//        const std::string &kernel_ = "gaussian",
//        int step = 1,
//        int verb = 0);
        
int main() {
    fmat I;
    fmat Is;

    I.load("lena.txt", raw_ascii);
    cout << "image shape: " << I.n_rows << " "<< I.n_cols << endl;
    
    vl_imsmooth(I, Is, 5.0);
    cout << "smooth image shape: " << Is.n_rows << " "<< Is.n_cols << endl;

    fmat Is2 = Is.submat(0,0,9,9);
    Is2.save("lena_smooth_cpp.txt", raw_ascii);
    return 0;
}
