#ifndef VLFEAT_HPP_
#define VLFEAT_HPP_

#include <armadillo>
#include <cstdlib>
#include <iostream>
#include <string>

inline void m_assert(bool condition, const std::string &message="") {
    if (!condition) {
        std::cout << message << std::endl;
        std::exit(0);
    }
}

void vl_sift(const arma::fmat &I, arma::mat &f, arma::fmat &d,
        int octaves = -1,
        int levels = -1,
        int first_octave = -1,
        double peak_thresh = -1,
        double edge_thresh = -1,
        double norm_thresh = -1,
        double magnif = -1,
        double window_size = -1,
        bool orientations = false,
        bool float_descriptors = false,
        int verbose = 0);


#endif /* VLFEAT_HPP_ */
