#ifndef THEBLOCK_H
#define THEBLOCK_H

#include "Hamiltonian.h"

class TheBlock;
class FinalSuperblock;

struct stepData
{
    Hamiltonian ham;                             // model Hamiltonian paramters
    bool sweepingEast,                 // currently sweeping from west to east?
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
        effectiveHams blockParts; // stores effective Hamiltonians for this site

        TheBlock() {};
        TheBlock(const Hamiltonian& ham, bool westSide);
        TheBlock(const effectiveHams& blockParts);
        TheBlock nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                           double& cumulativeTruncationError,
                           TheBlock* nextCompBlock = NULL);
                                                     // performs each DMRG step
        rmMatrixX_t projectBlock(const rmMatrixX_t& blockOp),
                          // projects system block into truncated DM eigenbasis
                    projectBASCoupling(const rmMatrixX_t& blockOp,
                                      const rmMatrixX_t& siteOp),
     // projects term coupling system block to left-hand ("adjacent") free site
                    projectFreeSite(const MatrixD_t& freeSiteOp);
                                                     // projects free site term
        FinalSuperblock createHSuperFinal(const stepData& data,
                                          rmMatrixX_t& psiGround, int skips)
                                          const;
               // identical to member functions above except for variable types
        #ifdef differentScalars
        obsMatrixX_t obsProjectBlock(const obsMatrixX_t& sysOp),
                     obsProjectBASCoupling(const obsMatrixX_t& blockOp,
                                          const obsMatrixX_t& siteOp),
                     obsProjectFreeSite(const obsMatrixD_t& lFreeSite);
        #endif
    
    private:
        rmMatrixX_t primeToRhoBasis;                  // change-of-basis matrix
        
        double lanczos(const stepData& data, rmMatrixX_t& seed) const;
     // changes input seed to ground eigenvector - make sure seed is normalized
        rmMatrixX_t createPrimeToRhoBasis(const MatrixX_t& rho, int mMax,
                                          double& cumulativeTruncationError)
                                          const,
           // create projection matrix to truncated Hilbert space of next block
                    createNewHS(const Hamiltonian& ham,
                                int freeSiteDistFromWestEnd);
                       // project expanded system block down to truncated basis
        void projectCouplingOperators(effectiveHams& newBlockParts,
                                      const Hamiltonian& ham);
                          // project coupling operators down to truncated basis
};

#endif
