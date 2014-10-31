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
        int firstOctave = -1,
        double peakThresh = -1,
        double edgeThresh = -1,
        double normThresh = -1,
        double magnif = -1,
        double windowSize = -1,
        bool orientations = false,
        bool floatDescriptors = false,
        int verbose = 0);


#endif /* VLFEAT_HPP_ */
