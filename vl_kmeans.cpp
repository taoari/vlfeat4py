/** @file   vl_kmeans.c
 ** @brief  vl_kmeans MEX definition.
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include <vl/kmeans.h>
//#include <mexutils.h>
#include "vlfeat.hpp"
#include <string.h>
#include <stdio.h>

void vl_kmeans(const arma::fmat &X, arma::fmat &Y, // for assignments and energy
        int numCenters,
        const std::string &algorithm_,
        const std::string &distance_,
        const std::string &initialization_,
        double minEnergyVariation,
        int numRepetitions,
        int numTrees,
        int maxNumComparisons,
        int maxNumIterations,
        int verbose)
{
    vl_size dimension ;
    vl_size numData ;
    vl_type dataType = VL_TYPE_FLOAT ; // TODO: only suport float

    void const * data = NULL ;

    VlKMeansAlgorithm algorithm = VlKMeansLloyd ;
    VlVectorComparisonType distance = VlDistanceL2 ;
//    vl_size maxNumIterations = 100 ;
//    vl_size numRepetitions = 1 ;
//    double minEnergyVariation = -1 ;
    double energy ;
    int verbosity = 0 ;
    VlKMeansInitialization initialization;
//    vl_size maxNumComparisons = 100 ;
//    vl_size numTrees = 3;

    VlKMeans * kmeans ;
    
    /* Input */
    // X : data of matrix DxN
    dimension = X.n_rows;
    numData = X.n_cols;
    // numCenters
    m_assert(numCenters >= 1 && numCenters <= numData,
            "NUMCENTERS must be a positive integer not greater "
            "than the number of data.");
    /* Options */
    // algorithm_ -> algorithm
    m_assert(algorithm_ == "lloyd" || algorithm_ == "elkan" || algorithm_ == "ann",
            "ALGORITHM must be one of LLOYD, ELKAN, ANN");
    if (algorithm_ == "lloyd")
        algorithm = VlKMeansLloyd;
    else if (algorithm_ == "elkan")
        algorithm = VlKMeansElkan;
    else if (algorithm_ == "ann")
        algorithm = VlKMeansANN;
    // initialization_ -> initialization
    m_assert(initialization_ == "plusplus" || initialization_ == "++" || initialization_ == "randsel",
            "INITIALIZATION must be one of PLUSPLUS/++, RANDSEL");
    if (initialization_ == "plusplus" || initialization_ == "++")
        initialization = VlKMeansPlusPlus;
    else if (initialization_ == "randsel")
        initialization = VlKMeansRandomSelection;
    // distance_ -> distance
    m_assert(distance_ == "l2" || distance_ == "l1" || distance_ == "chi2",
            "DISTANCE must be one of L2, L1, CHI2");
    if (distance_ == "l2")
        distance = VlDistanceL2;
    else if (distance_ == "l1")
        distance = VlDistanceL1;
    else if (distance_ == "chi2")
        distance = VlDistanceChi2;
    // numRepetitions
    m_assert(numRepetitions >= 1, 
            "NUMREPETITIONS must be larger than or equal to 1.");
    // numTrees
    m_assert(numTrees >= 1, 
            "NUMTREES must be larger than or equal to 1.");
    // maxNumComparisons
    m_assert(maxNumComparisons >= 0,
            "NUMCOMPARISONS must be larger than or equal to 0.");
    // maxNumIterations
    m_assert(maxNumIterations > 0, 
            "MAXNUMITERATIONS must be a non-negative integer");
    // minEnergyVariation
    m_assert(minEnergyVariation == -1 || minEnergyVariation > 0,
            "MINENERGYVARIATION must be a non-negative");
    // verbose -> verbosity
    verbosity = verbose;



//enum {
//  opt_max_num_iterations,
//  opt_algorithm,
//  opt_distance,
//  opt_initialization,
//  opt_num_repetitions,
//  opt_verbose,
//  opt_num_comparisons,
//  opt_min_energy_variation,
//  opt_num_trees,
//  opt_multithreading
//} ;

//enum {
//  INIT_RANDSEL,
//  INIT_PLUSPLUS
//} ;

//vlmxOption  options [] = {
//  {"MaxNumIterations",  1,   opt_max_num_iterations  },
//  {"Algorithm",         1,   opt_algorithm           },
//  {"Distance",          1,   opt_distance            },
//  {"Verbose",           0,   opt_verbose             },
//  {"NumRepetitions",    1,   opt_num_repetitions,    },
//  {"Initialization",    1,   opt_initialization      },
//  {"Initialisation",    1,   opt_initialization      }, /* UK spelling */
//  {"NumTrees",          1,   opt_num_trees           },
//  {"MaxNumComparisons", 1,   opt_num_comparisons     },
//  {"MinEnergyVariation",1,   opt_min_energy_variation},
//  {0,                   0,   0                       }
//} ;

///* driver */
//void
//mexFunction (int nout, mxArray * out[], int nin, const mxArray * in[])
//{

//  enum {IN_DATA = 0, IN_NUMCENTERS, IN_END} ;
//  enum {OUT_CENTERS = 0, OUT_ASSIGNMENTS, OUT_ENERGY} ;

//  int opt ;
//  int next = IN_END ;
//  mxArray const  *optarg ;

//  vl_size numCenters ;
//  vl_size dimension ;
//  vl_size numData ;

//  void const * data = NULL ;

//  VlKMeansAlgorithm algorithm = VlKMeansLloyd ;
//  VlVectorComparisonType distance = VlDistanceL2 ;
//  vl_size maxNumIterations = 100 ;
//  vl_size numRepetitions = 1 ;
//  double minEnergyVariation = -1 ;
//  double energy ;
//  int verbosity = 0 ;
//  int initialization = INIT_PLUSPLUS ;
//  vl_size maxNumComparisons = 100 ;
//  vl_size numTrees = 3;

//  vl_type dataType ;
//  mxClassID classID ;

//  VlKMeans * kmeans ;

//  VL_USE_MATLAB_ENV ;

//  /* -----------------------------------------------------------------
//   *                                               Check the arguments
//   * -------------------------------------------------------------- */

//  if (nin < 2) {
//    vlmxError (vlmxErrInvalidArgument,
//              "At least two arguments required.");
//  }
//  else if (nout > 3) {
//    vlmxError (vlmxErrInvalidArgument,
//              "Too many output arguments.");
//  }

//  classID = mxGetClassID (IN(DATA)) ;
//  switch (classID) {
//    case mxSINGLE_CLASS: dataType = VL_TYPE_FLOAT ; break ;
//    case mxDOUBLE_CLASS: dataType = VL_TYPE_DOUBLE ; break ;
//    default:
//      vlmxError (vlmxErrInvalidArgument,
//                "DATA must be of class SINGLE or DOUBLE") ;
//      abort() ;
//  }

//  dimension = mxGetM (IN(DATA)) ;
//  numData = mxGetN (IN(DATA)) ;

//  if (dimension == 0) {
//    vlmxError (vlmxErrInvalidArgument, "SIZE(DATA,1) is zero") ;
//  }

//  if (!vlmxIsPlainScalar(IN(NUMCENTERS)) ||
//      (numCenters = (vl_size) mxGetScalar(IN(NUMCENTERS))) < 1  ||
//      numCenters > numData) {
//    vlmxError (vlmxErrInvalidArgument,
//              "NUMCENTERS must be a positive integer not greater "
//              "than the number of data.") ;
//  }

//  while ((opt = vlmxNextOption (in, nin, options, &next, &optarg)) >= 0) {
//    char buf [1024] ;

//    switch (opt) {

//      case opt_verbose :
//        ++ verbosity ;
//        break ;

//      case opt_max_num_iterations :
//        if (!vlmxIsPlainScalar(optarg) || mxGetScalar(optarg) < 0) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "MAXNUMITERATIONS must be a non-negative integer scalar") ;
//        }
//        maxNumIterations = (vl_size) mxGetScalar(optarg) ;
//        break ;
//        
//      case opt_min_energy_variation :
//        if (!vlmxIsPlainScalar(optarg) || mxGetScalar(optarg) < 0) {
//          vlmxError (vlmxErrInvalidArgument,
//                     "MINENERGYVARIATION must be a non-negative scalar") ;
//        }
//        minEnergyVariation = mxGetScalar(optarg) ;
//        break ;

//      case opt_algorithm :
//        if (!vlmxIsString (optarg, -1)) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "ALGORITHM must be a string.") ;
//        }
//        if (mxGetString (optarg, buf, sizeof(buf))) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "ALGORITHM argument too long.") ;
//        }
//        if (vlmxCompareStringsI("lloyd", buf) == 0) {
//          algorithm = VlKMeansLloyd ;
//        } else if (vlmxCompareStringsI("elkan", buf) == 0) {
//          algorithm = VlKMeansElkan ;
//        } else if (vlmxCompareStringsI("ann", buf) == 0) {
//          algorithm = VlKMeansANN ;
//        } else {
//          vlmxError (vlmxErrInvalidArgument,
//                    "Invalid value %s for ALGORITHM", buf) ;
//        }
//        break ;

//      case opt_initialization :
//        if (!vlmxIsString (optarg, -1)) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "INITLAIZATION must be a string.") ;
//        }
//        if (mxGetString (optarg, buf, sizeof(buf))) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "INITIALIZATION argument too long.") ;
//        }
//        if (vlmxCompareStringsI("plusplus", buf) == 0 ||
//            vlmxCompareStringsI("++", buf) == 0) {
//          initialization = VlKMeansPlusPlus ;
//        } else if (vlmxCompareStringsI("randsel", buf) == 0) {
//          initialization = VlKMeansRandomSelection ;
//        } else {
//          vlmxError (vlmxErrInvalidArgument,
//                    "Invalid value %s for INITIALISATION.", buf) ;
//        }
//        break ;

//      case opt_distance :
//        if (!vlmxIsString (optarg, -1)) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "DISTANCE must be a string.") ;
//        }
//        if (mxGetString (optarg, buf, sizeof(buf))) {
//          vlmxError (vlmxErrInvalidArgument,
//                    "DISTANCE argument too long.") ;
//        }
//        if (vlmxCompareStringsI("l2", buf) == 0) {
//          distance = VlDistanceL2 ;
//        } else if (vlmxCompareStringsI("l1", buf) == 0) {
//          distance = VlDistanceL1 ;
//        } else if (vlmxCompareStringsI("chi2", buf) == 0) {
//          distance = VlDistanceChi2 ;
//        } else {
//          vlmxError (vlmxErrInvalidArgument,
//                    "Invalid value %s for DISTANCE", buf) ;
//        }
//        break ;

//      case opt_num_repetitions :
//        if (!vlmxIsPlainScalar (optarg)) {
//          vlmxError (vlmxErrInvalidArgument,
//                     "NUMREPETITIONS must be a scalar.") ;
//        }
//        if (mxGetScalar (optarg) < 1) {
//          vlmxError (vlmxErrInvalidArgument,
//                     "NUMREPETITIONS must be larger than or equal to 1.") ;
//        }
//        numRepetitions = (vl_size) mxGetScalar (optarg) ;
//        break ;

//       case opt_num_trees :
//            if (!vlmxIsPlainScalar (optarg)) {
//              vlmxError (vlmxErrInvalidArgument,
//                     "NUMTREES must be a scalar.") ;
//            }
//            if (mxGetScalar (optarg) < 1) {
//              vlmxError (vlmxErrInvalidArgument,
//                    "NUMTREES must be larger than or equal to 1.") ;
//            }
//            numTrees = (vl_size) mxGetScalar (optarg) ;
//         break;


//       case opt_num_comparisons :
//            if (!vlmxIsPlainScalar (optarg)) {
//              vlmxError (vlmxErrInvalidArgument,
//                     "NUMCOMPARISONS must be a scalar.") ;
//            }
//            if (mxGetScalar (optarg) < 0) {
//              vlmxError (vlmxErrInvalidArgument,
//                    "NUMCOMPARISONS must be larger than or equal to 0.") ;
//            }
//            maxNumComparisons = (vl_size) mxGetScalar (optarg) ;
//         break;

//      default :
//        abort() ;
//        break ;
//    }
//  }

  /* -----------------------------------------------------------------
   *                                                        Do the job
   * -------------------------------------------------------------- */

  data = X.memptr(); // mxGetPr(IN(DATA)) ;

  kmeans = vl_kmeans_new (dataType, distance) ;

  vl_kmeans_set_verbosity (kmeans, verbosity) ;
  vl_kmeans_set_num_repetitions (kmeans, numRepetitions) ;
  vl_kmeans_set_algorithm (kmeans, algorithm) ;
  vl_kmeans_set_initialization (kmeans, initialization) ;
  vl_kmeans_set_max_num_iterations (kmeans, maxNumIterations) ;
  vl_kmeans_set_max_num_comparisons (kmeans, maxNumComparisons) ;
  vl_kmeans_set_num_trees (kmeans, numTrees);
  
  if (minEnergyVariation >= 0) {
    printf("%f\n\n\n",minEnergyVariation);
    vl_kmeans_set_min_energy_variation (kmeans, minEnergyVariation) ;
  }

  if (verbosity) {
    char const * algorithmName = 0 ;
    char const * initializationName = 0 ;

    switch (vl_kmeans_get_algorithm(kmeans)) {
      case VlKMeansLloyd: algorithmName = "Lloyd" ; break ;
      case VlKMeansElkan: algorithmName = "Elkan" ; break ;
      case VlKMeansANN:   algorithmName = "ANN" ; break ;
      default : abort() ;
    }
    switch (vl_kmeans_get_initialization(kmeans)) {
      case VlKMeansPlusPlus : initializationName = "plusplus" ; break ;
      case VlKMeansRandomSelection : initializationName = "randsel" ; break ;
      default: abort() ;
    }
    printf("kmeans: Initialization = %s\n", initializationName) ;
    printf("kmeans: Algorithm = %s\n", algorithmName) ;
    printf("kmeans: MaxNumIterations = %d\n", (int)vl_kmeans_get_max_num_iterations(kmeans)) ;
    printf("kmeans: MinEnergyVariation = %f\n", vl_kmeans_get_min_energy_variation(kmeans)) ;
    printf("kmeans: NumRepetitions = %d\n", (int)vl_kmeans_get_num_repetitions(kmeans)) ;
    printf("kmeans: data type = %s\n", vl_get_type_name(vl_kmeans_get_data_type(kmeans))) ;
    printf("kmeans: distance = %s\n", vl_get_vector_comparison_type_name(vl_kmeans_get_distance(kmeans))) ;
    printf("kmeans: data dimension = %d\n", (int)dimension) ;
    printf("kmeans: num. data points = %d\n", (int)numData) ;
    printf("kmeans: num. centers = %d\n", numCenters) ;
    printf("kmeans: max num. comparisons = %d\n", maxNumComparisons) ;
    printf("kmeans: num. trees = %d\n", numTrees) ;
    printf("\n") ;
  }

  /* -------------------------------------------------------------- */
  /*                                    Clustering and quantization */
  /* -------------------------------------------------------------- */

    // TODO: ! cause memory leak
  // energy = vl_kmeans_cluster(kmeans, data, dimension, numData, numCenters) ;

  /* copy centers */
  // OUT(CENTERS) = mxCreateNumericMatrix (dimension, numCenters, classID, mxREAL) ;
  Y = arma::fmat(dimension, numCenters);
  memcpy (Y.memptr(), // mxGetData(OUT(CENTERS)),
          vl_kmeans_get_centers (kmeans),
          vl_get_type_size (dataType) * dimension * vl_kmeans_get_num_centers(kmeans)) ;

// TODO: return assignments and energy
//  /* optionally qunatize */
//  if (nout > 1) {
//    vl_uindex j ;
//    vl_uint32 * assignments  ;
//    OUT(ASSIGNMENTS) = mxCreateNumericMatrix (1, numData, mxUINT32_CLASS, mxREAL) ;
//    assignments = mxGetData (OUT(ASSIGNMENTS)) ;

//    vl_kmeans_quantize (kmeans, assignments, NULL, data, numData) ;

//    /* use MATLAB indexing convention */
//    for (j = 0 ; j < numData ; ++j) { assignments[j] += 1 ; }
//  }

//  /* optionally return energy */
//  if (nout > 2) {
//    OUT(ENERGY) = vlmxCreatePlainScalar (energy) ;
//  }

  vl_kmeans_delete (kmeans) ;
}
