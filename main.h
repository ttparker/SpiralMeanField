#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <unsupported/Eigen/KroneckerProduct>
#include "GlobalHamiltonianParameters.h"

#ifdef realHamiltonian
    typedef double hamScalar;
#elif defined(complexHamiltonian)
    typedef std::complex<double> hamScalar;
#endif

#if defined(realHamiltonian) && defined(realObservables)
    #define obsRe
    typedef double obsScalar;
#elif defined(complexObservables) || defined(complexHamiltonian)
    #define obsRe std::real
    typedef std::complex<double> obsScalar;
#endif

#define obsId_d Matrix<obsScalar, d, d>::Identity()
#define obsId(size) Matrix<obsScalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(size, size)

typedef Eigen::Matrix<hamScalar, d, d> MatrixD_t;
typedef Eigen::Matrix<hamScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixX_t;
typedef Eigen::Matrix<hamScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    rmMatrixX_t;
typedef Eigen::Matrix<hamScalar, Eigen::Dynamic, 1> VectorX_t;
typedef Eigen::Matrix<obsScalar, d, d> obsMatrixD_t;
typedef Eigen::Matrix<obsScalar, Eigen::Dynamic, Eigen::Dynamic> obsMatrixX_t;

#endif
