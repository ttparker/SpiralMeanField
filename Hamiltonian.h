#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"

struct effectiveHams       // each site's effective block and bond Hamiltonians
{
    int l,                // site at the end of the block (i.e. block size - 1)
        m;                                  // number of states stored in block
    rmMatrixX_t hS;                                        // block Hamiltonian
    std::vector<MatrixX_t> off0RhoBasisH2,
                           off1RhoBasisH2;
            // density-matrix-basis coupling operators - "off" means the offset
            // between this block, in which the operator is represented, and
            // the site on which it acts
};

class TheBlock;

class Hamiltonian
{
    public:
        Hamiltonian();
        void setParams(const std::vector<double>& couplingConstants, int lSys,
                       double k);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        int lSys;                                      // current system length
        std::vector<double> couplingConstants;
        std::vector<MatrixD_t, Eigen::aligned_allocator<MatrixD_t>> siteBasisH2;
                                     // site-basis coupling operators - the
                                     // independent ones should be listed first
        double h;                           // external magnetic field strength
        std::vector<double> cosList,
                            sinList;
        MatrixD_t sigmax;                             // for the external field
        
        rmMatrixX_t EDBASCoupling(int j,
                                  const std::vector<MatrixX_t>& rhoBasisH2) const;
                  // exact-diagonalization stage block-adjacent-site coupling -
                  // j gives the (j-1)th coupling constant
        MatrixD_t h1(int sitesFromWestEnd) const;
                        // returns the external magnetic field as a function of
                        // position, measured from the initially leftmost site
        Eigen::Matrix<hamScalar, Eigen::Dynamic, 1>
            act(const effectiveHams& blockParts,
                const effectiveHams& compBlockParts, rmMatrixX_t x,
                bool sweepingEast) const;
        rmMatrixX_t sysBlockRSiteCoupling(const MatrixX_t& hSBond,
                                          const MatrixD_t& hSiteBond,
                                          const rmMatrixX_t& x, int mSIn,
                                          int mEIn) const,
                    sysBlockLSiteCoupling(const MatrixX_t& hSBond,
                                          const MatrixD_t& hSiteBond,
                                          const rmMatrixX_t& cycledx, int mSIn,
                                          int MEIn) const;
        Eigen::Matrix<hamScalar, 1, Eigen::Dynamic>
            freeSiteCouplingTerm(const MatrixD_t& lFreeSiteBond,
                                 const MatrixD_t& rFreeSiteBond,
                                 const rmMatrixX_t& RHS, int mEIn) const;
        rmMatrixX_t envBlockLSiteCoupling(const MatrixX_t& hEBond,
                                          const MatrixD_t& hSiteBond,
                                          const rmMatrixX_t& cycledx, int mSIn,
                                          int mEIn) const,
                    envBlockRSiteCoupling(const MatrixX_t& hEBond,
                                          const MatrixD_t& hSiteBond,
                                          const rmMatrixX_t& cycledx, int mSIn,
                                          int MEIn) const,
                    projectedBASCouplings(TheBlock * block) const;
                                    // specify the block-adjacent-free-site
                                    // couplings to be projected into new basis
    
    friend class TheBlock;
};

#endif
