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

void vl_imsmooth(const arma::fmat &I, arma::fmat &Is, 
        double sigma,
        const std::string &padding_ = "continuity",
        const std::string &kernel_ = "gaussian",
        int subsample = 1,
        int verbose = 0);

void vl_imsmooth(const arma::mat &I, arma::mat &Is, 
        double sigma,
        const std::string &padding_ = "continuity",
        const std::string &kernel_ = "gaussian",
        int subsample = 1,
        int verbose = 0);

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

void vl_dsift(const arma::fmat &I, arma::mat &f, arma::fmat &d,
        const arma::mat &bounds_ = arma::mat(),
        const arma::mat &step_ = arma::mat(),
        const arma::mat &size_ = arma::mat(),
        const arma::mat &geometry_ = arma::mat(),
        bool fast = true,
        bool norm = true,
        double windowSize = -1,
        bool floatDescriptors = false,
        int verbose = 0);

void vl_homkermap(const arma::fmat &XX, arma::fmat &VV, int n,
        const std::string & kernel = "kchi2",
        const std::string & window = "rectangular",
        double gamma = 1.0,
        double period = -1);

#endif /* VLFEAT_HPP_ */
