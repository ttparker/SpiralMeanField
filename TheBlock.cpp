#include "FinalSuperblock.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

TheBlock::TheBlock(const Hamiltonian& ham, bool westSide)
{
    blockParts.l = 0;
    blockParts.m = d;
    blockParts.hS = (westSide ? ham.h1(0) : ham.h1(ham.lSys - 1));
    blockParts.off0RhoBasisH2
        .assign(ham.siteBasisH2.begin(),
                ham.siteBasisH2.begin() + nIndepCouplingOperators);
};

TheBlock::TheBlock(const effectiveHams& blockParts) : blockParts(blockParts) {};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                             double& cumulativeTruncationError,
                             TheBlock* nextCompBlock)
{
    int md = blockParts.m * d;
    effectiveHams newBlockParts;                  // data for next system block
    newBlockParts.l = blockParts.l + 1;
    if(data.exactDiag)
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    {
        newBlockParts.m = md;
        newBlockParts.off0RhoBasisH2.reserve(nIndepCouplingOperators);
        newBlockParts.off1RhoBasisH2.reserve(nIndepCouplingOperators);
        for(int i = 0; i < nIndepCouplingOperators; i++)
        {
            newBlockParts.off0RhoBasisH2
                .push_back(kp(MatrixXd::Identity(blockParts.m, blockParts.m),
                              data.ham.siteBasisH2[i]));
            newBlockParts.off1RhoBasisH2
                .push_back(kp(blockParts.off0RhoBasisH2[i],
                              MatrixD_t::Identity()));
        };
        effectiveHams newCompBlockParts = newBlockParts;
                                             // data for next environment block
        newBlockParts.hS.noalias()
            =  kp(blockParts.hS, Matrix<double, d, d>::Identity())
             + data.ham.EDBASCoupling(1, blockParts.off0RhoBasisH2)
             + kp(MatrixXd::Identity(blockParts.m, blockParts.m),
                  data.ham.h1(blockParts.l + 1));
        newCompBlockParts.hS.noalias()
            =  kp(data.compBlock -> blockParts.hS,
                  Matrix<double, d, d>::Identity())
             + data.ham.EDBASCoupling(1, data.compBlock
                                         -> blockParts.off0RhoBasisH2)
             + kp(MatrixXd::Identity(data.compBlock -> blockParts.m,
                                     data.compBlock -> blockParts.m),
                  data.ham.h1(data.ham.lSys - 2 - data.compBlock -> blockParts.l));
        if(blockParts.l != 0)
        {
            newBlockParts.hS.noalias()
                += data.ham.EDBASCoupling(2, blockParts.off1RhoBasisH2);
            newCompBlockParts.hS.noalias()
                += data.ham.EDBASCoupling(2, data.compBlock
                                             -> blockParts.off1RhoBasisH2);
        };
        *nextCompBlock = TheBlock(newCompBlockParts);
        return TheBlock(newBlockParts);
    };
    lanczos(data, psiGround);                         // calculate ground state
    int compm = data.compBlock -> blockParts.m;
    psiGround.resize(md, compm * d);
    
    // project the expanded system block into a new block:
    primeToRhoBasis.noalias()
        = createPrimeToRhoBasis(psiGround * psiGround.adjoint(), data.mMax,
                                cumulativeTruncationError);
    newBlockParts.m = primeToRhoBasis.cols();
                        // number of states kept in next truncated system block
    int lFreeSiteDistFromWestEnd,
        rFreeSiteDistFromWestEnd;
    if(data.sweepingEast)
    {
        lFreeSiteDistFromWestEnd = blockParts.l + 1;
        rFreeSiteDistFromWestEnd = data.ham.lSys - 2
                                   - data.compBlock -> blockParts.l;
    }
    else
    {
        lFreeSiteDistFromWestEnd = data.ham.lSys - 2 - blockParts.l;
        rFreeSiteDistFromWestEnd = data.compBlock -> blockParts.l + 1;
    };
    newBlockParts.hS.noalias() = createNewHS(data.ham, lFreeSiteDistFromWestEnd);
    projectCouplingOperators(newBlockParts, data.ham);
    
    if(data.infiniteStage)
    // project the expanded environment block into a new block
    {
        if(nextCompBlock)           // excludes last odd-system-size iDMRG step
        {
            data.compBlock -> primeToRhoBasis.noalias()
                = createPrimeToRhoBasis(psiGround.adjoint() * psiGround,
                                        data.mMax, cumulativeTruncationError);
                                // construct environment change-of-basis matrix
            effectiveHams newCompBlockParts;  // data for new environment block
            newCompBlockParts.l = data.compBlock -> blockParts.l + 1;
            newCompBlockParts.m = data.compBlock -> primeToRhoBasis.cols();
            newCompBlockParts.hS.noalias()
                = data.compBlock
                    -> createNewHS(data.ham, rFreeSiteDistFromWestEnd);
            data.compBlock
                -> projectCouplingOperators(newCompBlockParts, data.ham);
            *nextCompBlock = TheBlock(newCompBlockParts);
        };
    }
    else                   // modify psiGround to predict the next ground state
    {
        for(int sPrimeIndex = 0; sPrimeIndex < md; sPrimeIndex++)
                // transpose the environment block and the right-hand free site
        {
            rmMatrixX_t ePrime = psiGround.row(sPrimeIndex);
            ePrime.resize(compm, d);
            ePrime.transposeInPlace();
            ePrime.resize(1, d * compm);
            psiGround.row(sPrimeIndex) = ePrime;
        };
        psiGround = primeToRhoBasis.adjoint() * psiGround; 
                                      // change the expanded system block basis
        psiGround.resize(newBlockParts.m * d, compm);
        psiGround *= data.beforeCompBlock -> primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(newBlockParts.m * d
                         * data.beforeCompBlock -> blockParts.m * d, 1);
    };
    return TheBlock(newBlockParts); // save expanded-block operators in new basis
};

rmMatrixX_t TheBlock::createPrimeToRhoBasis(const MatrixX_t& rho, int mMax,
                                            double& cumulativeTruncationError)
                                            const
{
    SelfAdjointEigenSolver<MatrixX_t> rhoSolver(rho);
                                             // find density matrix eigenstates
    int md = blockParts.m * d,
        evecsToKeep;
    if(md <= mMax)
        evecsToKeep = md;
    else
    {
        int firstKeptEval = md - mMax;
        for(; firstKeptEval < md
              && (rhoSolver.eigenvalues()(firstKeptEval) == 0
                  ||   (  rhoSolver.eigenvalues()(firstKeptEval)
                        - rhoSolver.eigenvalues()(firstKeptEval - 1))
                       / std::abs(rhoSolver.eigenvalues()(firstKeptEval))
                     < degenerateDMCutoff);
              firstKeptEval++);
              // find the the max number of eigenvectors to keep that do not
              // terminate inside a degenerate eigenspace of the density matrix
        evecsToKeep = md - firstKeptEval;
        if(evecsToKeep == 0)
        {
            std::cerr << "More than mMax highest-weighted density-matrix "
                      << "eigenvectors are degenerate." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if(evecsToKeep != mMax)
            std::cout << "Warning: mMax truncation ends in a degenerate DM "
                      << "eigenspace, lowering cutoff to " << evecsToKeep
                      << " states." << std::endl;
    };
    cumulativeTruncationError
        += rhoSolver.eigenvalues().head(md - evecsToKeep).sum();
    return rhoSolver.eigenvectors().rightCols(evecsToKeep);
                                     // construct system change-of-basis matrix
};

rmMatrixX_t TheBlock::createNewHS(const Hamiltonian& ham,
                                  int freeSiteDistFromWestEnd)
{
    return  projectBlock(blockParts.hS)                    // system block term
          + ham.projectedBASCouplings(this)                    // coupling term
          + projectFreeSite(ham.h1(freeSiteDistFromWestEnd)); // free site term
};

rmMatrixX_t TheBlock::projectBlock(const rmMatrixX_t& blockOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    rmMatrixX_t oSO = blockOp * primeToRhoBasis;
    oSO.resize(blockParts.m * d, nextSiteM);
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    return primeToRhoBasis.adjoint() * oSO;
};

rmMatrixX_t TheBlock::projectBASCoupling(const rmMatrixX_t& blockOp,
                                         const rmMatrixX_t& siteOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    rmMatrixX_t oDagHSite = oDag * siteOp,
                hSysO = blockOp * primeToRhoBasis;
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    oDagHSite.resize(nextSiteM, blockParts.m * d);
    hSysO.resize(blockParts.m * d, nextSiteM);
    return oDagHSite * hSysO;
};

rmMatrixX_t TheBlock::projectFreeSite(const MatrixD_t& freeSiteOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    rmMatrixX_t oDagH = oDag * freeSiteOp;
    oDagH.resize(nextSiteM, blockParts.m * d);
    return oDagH * primeToRhoBasis;
};

void TheBlock::projectCouplingOperators(effectiveHams& newBlockParts,
                                        const Hamiltonian& ham)
{
    newBlockParts.off0RhoBasisH2.reserve(nIndepCouplingOperators);
    newBlockParts.off1RhoBasisH2.reserve(nIndepCouplingOperators);
    for(int i = 0; i < nIndepCouplingOperators; i++)
                     // project coupling operators into new DM eigenspace basis
    {
        newBlockParts.off0RhoBasisH2
            .push_back(projectFreeSite(ham.siteBasisH2[i]));
        newBlockParts.off1RhoBasisH2
            .push_back(projectBlock(blockParts.off0RhoBasisH2[i]));
    };
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            rmMatrixX_t& psiGround, int skips)
                                            const
{
    double gsEnergy = lanczos(data, psiGround);       // calculate ground state
    return FinalSuperblock(gsEnergy, data.ham.lSys, psiGround, blockParts.m,
                           data.compBlock -> blockParts.m, skips);
};

#ifdef differentScalars
obsMatrixX_t TheBlock::obsProjectBlock(const obsMatrixX_t& sysOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    obsMatrixX_t oSO = sysOp * primeToRhoBasis;
    oSO.resize(blockParts.m * d, nextSiteM);
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    return primeToRhoBasis.adjoint() * oSO;
};

obsMatrixX_t TheBlock::obsProjectNNCoupling(const obsMatrixX_t& blockOp,
                                            const obsMatrixX_t& siteOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    obsMatrixX_t oDagHSite = oDag * siteOp,
                 hSysO = blockOp * primeToRhoBasis;
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    oDagHSite.resize(nextSiteM, blockParts.m * d);
    hSysO.resize(blockParts.m * d, nextSiteM);
    return oDagHSite * hSysO;
};

obsMatrixX_t TheBlock::obsProjectFreeSite(const obsMatrixD_t& lFreeSite)
{
    int nextSiteM = primeToRhoBasis.cols();
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    obsMatrixX_t oDagH1 = oDag * lFreeSite;
    oDagH1.resize(nextSiteM, blockParts.m * d);
    return oDagH1 * primeToRhoBasis;
};
#endif
