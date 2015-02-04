#include "TheBlock.h"
#include "GlobalPrecisionParameters.h"

#ifdef realHamiltonian
    #define re
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

double TheBlock::lanczos(const stepData& data, rmMatrixX_t& seed) const
{
    int hamDimension = blockParts.m * d * data.compBlock -> blockParts.m * d;
    const int minIters = std::min(hamDimension, globalMinLancIters),
              maxIters = std::min(hamDimension, globalMaxLancIters);
    std::vector<double> a,
                        b;
    a.reserve(minIters);
    b.reserve(minIters);
    MatrixX_t basisVecs = seed;
    VectorX_t x = data.ham.act(blockParts, data.compBlock -> blockParts,
                               basisVecs, data.sweepingEast);
    
    
    double oldGSE,
           GSE = re(seed.col(0).dot(x));
    a.push_back(GSE);


    b.push_back(0.);
    int i = 0;                                             // iteration counter
//    char JOBZ = 'V',                                 // define dstemr arguments
//         RANGE = 'I';
    int N = 1;
    
    
    
//    std::vector<double> W;
    
    
    
/*    std::vector<double> D,
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
        INFO;*/
    double gSEDiff;
          // change in ground state energy across subsequent Lanczos iterations
    do
    {
        i++;
        oldGSE = GSE;
        
        // Lanczos stage 1: Lanczos iteration
        x -= a[i - 1] * basisVecs.col(i - 1);
        b.push_back(x.norm());
        basisVecs.conservativeResize(NoChange, i + 1);
        basisVecs.col(i) = x / b[i];
        x.noalias() =  data.ham.act(blockParts, data.compBlock -> blockParts,
                                    basisVecs.col(i), data.sweepingEast)
                     - b[i] * basisVecs.col(i - 1);
        a.push_back(re(basisVecs.col(i).dot(x)));
        
        // Lanczos stage 2: diagonalize tridiagonal matrix
        N++;
/*        D = a;
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
                &LWORK, IWORK.data(), &LIWORK, &INFO); // calculate ground state*/
        
        
        
        MatrixX_t triDiag = MatrixX_t::Zero(i + 1, i + 1);
        triDiag(0, 0) = a[0];
        for(int j = 1; j <= i; j++)
        {
            triDiag(j, j) = a[j];
            triDiag(j - 1, j) = triDiag(j, j - 1) = b[j];
        };
        SelfAdjointEigenSolver<MatrixX_t> solver(triDiag);
        GSE = solver.eigenvalues()(0);
        VectorX_t Z = solver.eigenvectors().col(0);
        seed = (basisVecs * Z).normalized();
        gSEDiff = std::abs((GSE - oldGSE) / oldGSE);
    } while(N < minIters || (N < maxIters && gSEDiff > data.lancTolerance));
    if(N == maxIters && gSEDiff > data.lancTolerance)
                          // check if last iteration converges to an eigenstate
    {
        std::cout << "Warning: final Lanczos iteration reached. The percent "
                  << "difference between the last two interations' ground state "
                  << "energies is " << gSEDiff << std::endl;
        if(gSEDiff > fallbackLancTolerance)
        {
            std::cerr << "Lanczos algorithm failed to converge after "
                      << maxIters << " iterations." << std::endl;
            exit(EXIT_FAILURE);
        };
    };
    
    std::cout << "Lanczos iterations: " << N << std::endl;
    
    
    return GSE;
};
