#include "TheBlock.h"
#include "GlobalPrecisionParameters.h"

#ifdef realHamiltonian
    #define re
    #include <cmath>
#endif
#ifdef complexHamiltonian
    #define re std::real
#endif

extern "C"
{
    void dstemr_(char* JOBZ, char* RANGE, int* N, double* D, double* E,
                 double* VL, double* VU, int* IL, int* IU, int* M, double* W,
                 double* Z, int* LDZ, int* NZC, int* ISUPPZ, bool* TRYRAC,
                 double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
};

using namespace Eigen;

double TheBlock::lanczos(const MatrixX_t& mat, rmMatrixX_t& seed,
                         double lancTolerance) const
{
    int matSize = mat.rows();
    if(matSize == 1)
        return re(mat(0, 0));
    const int minIters = std::min(matSize, globalMinLancIters),
              maxIters = std::min(matSize, globalMaxLancIters);
    std::vector<double> a,
                        b;
    a.reserve(minIters);
    b.reserve(minIters);
    MatrixX_t basisVecs = seed;
    VectorX_t x = mat * basisVecs;
    a.push_back(re(seed.col(0).dot(x)));
    b.push_back(0.);
    VectorX_t oldGS;
    int i = 0;                                             // iteration counter
    char JOBZ = 'V',                                 // define dstemr arguments
         RANGE = 'I';
    int N = 1;
    std::vector<double> D,
                        E;
    D.reserve(minIters);
    E.reserve(minIters);
    double VL,
           VU;
    int IL = 1,
        IU = 1,
        M;
    std::vector<double> W;
    W.reserve(minIters);
    VectorXd Z;
    int LDZ,
        NZC = 1,
        ISUPPZ[2];
    bool TRYRAC = true;
    double optLWORK;
    std::vector<double> WORK;
    int LWORK,
        optLIWORK;
    std::vector<int> IWORK;
    int LIWORK,
        INFO;
    double gStateDiff;
          // change in ground state vector across subsequent Lanczos iterations
    do
    {
        i++;
        oldGS = seed;
        
        // Lanczos stage 1: Lanczos iteration
        x -= a[i - 1] * basisVecs.col(i - 1);
        b.push_back(x.norm());
        basisVecs.conservativeResize(NoChange, i + 1);
        basisVecs.col(i) = x / b[i];
        x.noalias() = mat * basisVecs.col(i) - b[i] * basisVecs.col(i - 1);
        a.push_back(re(basisVecs.col(i).dot(x)));
        
        // Lanczos stage 2: diagonalize tridiagonal matrix
        N++;
        D = a;
        E.assign(b.begin() + 1, b.end());
        E.resize(N);
        W.resize(N);
        Z.resize(N);
        LDZ = N;
        LWORK = -1;
        LIWORK = -1;
        dstemr_(&JOBZ, &RANGE, &N, D.data(), E.data(), &VL, &VU, &IL, &IU, &M,
                W.data(), Z.data(), &LDZ, &NZC, ISUPPZ, &TRYRAC, &optLWORK,
                &LWORK, &optLIWORK, &LIWORK, &INFO);
                                     // query for optimal workspace allocations
        LWORK = int(optLWORK);
        WORK.resize(LWORK);
        LIWORK = optLIWORK;
        IWORK.resize(LIWORK);
        dstemr_(&JOBZ, &RANGE, &N, D.data(), E.data(), &VL, &VU, &IL, &IU, &M,
                W.data(), Z.data(), &LDZ, &NZC, ISUPPZ, &TRYRAC, WORK.data(),
                &LWORK, IWORK.data(), &LIWORK, &INFO); // calculate ground state
        seed = (basisVecs * Z).normalized();
        gStateDiff = std::abs(1 - std::abs(seed.col(0).dot(oldGS)));
    } while(N < minIters || (N < maxIters && gStateDiff > lancTolerance));
    if(N == maxIters && gStateDiff > lancTolerance)
                          // check if last iteration converges to an eigenstate
    {
        double gStateError = std::abs(1 - std::abs(seed.col(0)
                                                   .dot((mat * seed).col(0)
                                                        .normalized())));
        std::cout << "Warning: final Lanczos iteration reached. The inner "
                  << "product of the final approximate ground state and its "
                  << "normalized image differs from 1 by " << gStateError
                  << std::endl;
        if(gStateError > fallbackLancTolerance)
        {
            std::cerr << "Lanczos algorithm failed to converge after "
                      << maxIters << " iterations." << std::endl;
            exit(EXIT_FAILURE);
        };
    };
    return W.front();
};
