#ifndef THEBLOCK_H
#define THEBLOCK_H

#include "Hamiltonian.h"

typedef std::pair<MatrixX_t, MatrixX_t> matPair;
typedef std::pair<std::vector<MatrixX_t>, std::vector<MatrixX_t>> vecMatPair;

class FinalSuperblock;
class TheBlock;

struct stepData
{
    Hamiltonian ham;                             // model Hamiltonian paramters
    bool sweepingEast,                 // currently sweeping from east to west?
         exactDiag;             // close enough to edge to skip DMRG trucation?
    TheBlock* compBlock;         // complementary block on other side of system
    bool infiniteStage;           // currently performing infinite-system DMRG?
    double lancTolerance;  // max deviation from 1 of dot product of successive
                           // Lanczos iterations' ground state vectors
    int mMax;                              // max size of effective Hamiltonian
    TheBlock* beforeCompBlock;     // next smaller block than complementary one
};

class TheBlock
{
    public:
        int m;                              // number of states stored in block
        
        TheBlock(int m = 0, const MatrixX_t& hS = MatrixX_t(),
                 const vecMatPair& newRhoBasisH2s = vecMatPair(),
                 int l = 0);
        TheBlock(const Hamiltonian& ham, bool westSide);
        TheBlock nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                           double& cumulativeTruncationError,
                           TheBlock* nextCompBlock = NULL);
                                                     // performs each DMRG step
        FinalSuperblock createHSuperFinal(const stepData& data,
                                          rmMatrixX_t& psiGround, int skips)
                                          const;
        obsMatrixX_t obsChangeBasis(const obsMatrixX_t& mat) const;
                       // changes basis during calculation of observables stage
    
    private:
        MatrixX_t hS;                                      // block Hamiltonian
        MatrixX_t primeToRhoBasis;                    // change-of-basis matrix
        std::vector<MatrixX_t> off0RhoBasisH2,
                               off1RhoBasisH2;
            // density-matrix-basis coupling operators - "off" means the offset
            // between this block, in which the operator is represented, and
            // the site on which it acts
        int l;            // site at the end of the block (i.e. block size - 1)
        
        matPair createHprimes(const stepData& data) const;
        vecMatPair createNewRhoBasisH2s(const vecMatD_t& siteBasisH2,
                                     bool exactDiag, const TheBlock* block)
                                     const;
        double lanczos(const MatrixX_t& mat, rmMatrixX_t& seed,
                       double lancTolerance) const,
     // changes input seed to ground eigenvector - make sure seed is normalized
               solveHSuper(const matPair& hPrimes, const stepData& data,
                           rmMatrixX_t& psiGround) const;
        MatrixX_t createPrimeToRhoBasis(const MatrixX_t& rho, int mMax,
                                        double& cumulativeTruncationError),
           // create projection matrix to truncated Hilbert space of next block
                  changeBasis(const MatrixX_t& mat, const TheBlock* block) const;
                   // represents operators in the basis of the new system block
};

#endif
