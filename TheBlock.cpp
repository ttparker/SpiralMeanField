#include "FinalSuperblock.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixX_t& hS,
                   const vecMatPair& newRhoBasisH2s, int l)
    : m(m), hS(hS), off0RhoBasisH2(newRhoBasisH2s.first),
      off1RhoBasisH2(newRhoBasisH2s.second), l(l) {};

TheBlock::TheBlock(const Hamiltonian& ham, bool westSide)
    : m(d), hS(westSide ? ham.h1(0) : ham.h1(ham.lSys - 1)), l(0)
{
    off0RhoBasisH2.assign(ham.siteBasisH2.begin(),
                          ham.siteBasisH2.begin() + nIndepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                             double& cumulativeTruncationError,
                             TheBlock* nextCompBlock)
{
    matPair hPrimes = createHprimes(data);
                                      // expanded system and environment blocks
    int md = m * d;
    if(data.exactDiag)
    {
        *nextCompBlock = TheBlock(md, hPrimes.second,
                                  createNewRhoBasisH2s(data.ham.siteBasisH2,
                                                       true, this), l + 1);
        return TheBlock(md, hPrimes.first,
                        createNewRhoBasisH2s(data.ham.siteBasisH2, true, this),
                        l + 1);
    };
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    solveHSuper(hPrimes, data, psiGround);            // calculate ground state
    int compm = data.compBlock -> m;
    psiGround.resize(md, compm * d);
    primeToRhoBasis = createPrimeToRhoBasis(psiGround * psiGround.adjoint(),
                                            data.mMax,
                                            cumulativeTruncationError);
    int nextBlockM = primeToRhoBasis.cols();
                               // number of states kept in next truncated block
    if(data.infiniteStage)            // create corresponding environment block
    {
        data.compBlock -> primeToRhoBasis
            = createPrimeToRhoBasis(psiGround.adjoint() * psiGround, data.mMax,
                                    cumulativeTruncationError);
                                // construct environment change-of-basis matrix
        if(nextCompBlock)           // excludes last odd-system-size iDMRG step
            *nextCompBlock = TheBlock(data.compBlock -> primeToRhoBasis.cols(),
                                      changeBasis(hPrimes.second,
                                                  data.compBlock),
                                      createNewRhoBasisH2s(data.ham.siteBasisH2,
                                                           false,
                                                           data.compBlock),
                                      l + 1);
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
        psiGround.resize(nextBlockM * d, compm);
        psiGround *= data.beforeCompBlock -> primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(nextBlockM * d * data.beforeCompBlock -> m * d, 1);
    };
    return TheBlock(nextBlockM, changeBasis(hPrimes.first, this),
                    createNewRhoBasisH2s(data.ham.siteBasisH2, false, this),
                    l + 1);       // save expanded-block operators in new basis
};

matPair TheBlock::createHprimes(const stepData& data) const
{
    int hSprimeDistFromWestEnd,
        hEprimeDistFromWestEnd;
    if(data.sweepingEast)
    {
        hSprimeDistFromWestEnd = l + 1;
        hEprimeDistFromWestEnd = (data.infiniteStage ? data.ham.lSys - l - 2 :
                                  l + 2);
    }
    else
    {
        hSprimeDistFromWestEnd = data.ham.lSys - l - 2;
        hEprimeDistFromWestEnd = data.ham.lSys - l - 3;
    };
    MatrixX_t hSprime = kp(hS, Id_d)
                        + data.ham.blockAdjacentSiteJoin(1, off0RhoBasisH2)
                        + kp(Id(m), data.ham.h1(hSprimeDistFromWestEnd));
    if(l != 0)
        hSprime += data.ham.blockAdjacentSiteJoin(2, off1RhoBasisH2);
    MatrixX_t hEprime = kp(data.compBlock -> hS, Id_d)
                        + data.ham.blockAdjacentSiteJoin(1, data.compBlock
                                                            -> off0RhoBasisH2)
                        + kp(Id(data.compBlock -> m),
                             data.ham.h1(hEprimeDistFromWestEnd));
    if(data.compBlock -> l != 0)
        hEprime += data.ham.blockAdjacentSiteJoin(2, data.compBlock
                                                     -> off1RhoBasisH2);
    return std::make_pair(hSprime, hEprime);
};

vecMatPair TheBlock::createNewRhoBasisH2s(const vecMatD_t& siteBasisH2,
                                          bool exactDiag,
                                          const TheBlock* block) const
{
    std::vector<MatrixX_t> newOff0RhoBasisH2,
                           newOff1RhoBasisH2;
    newOff0RhoBasisH2.reserve(nIndepCouplingOperators);
    newOff1RhoBasisH2.reserve(nIndepCouplingOperators);
    for(int i = 0; i < nIndepCouplingOperators; i++)
    {
        newOff0RhoBasisH2.push_back(exactDiag ?
                                    kp(Id(m), siteBasisH2[i]) :
                                    changeBasis(kp(Id(m), siteBasisH2[i]),
                                                block));
        newOff1RhoBasisH2.push_back(exactDiag ?
                                    kp(off0RhoBasisH2[i], Id_d) :
                                    changeBasis(kp(block -> off0RhoBasisH2[i],
                                                   Id_d),
                                                block));
    };
    return std::make_pair(newOff0RhoBasisH2, newOff1RhoBasisH2);
};

double TheBlock::solveHSuper(const matPair& hPrimes, const stepData& data,
                             rmMatrixX_t& psiGround) const
{

    int compm = data.compBlock -> m;
    MatrixX_t hSuper = kp(hPrimes.first, Id(compm * d))
                       + data.ham.lBlockrSiteJoin(off0RhoBasisH2, compm)
                       + data.ham.siteSiteJoin(m, compm)
                       + data.ham.lSiterBlockJoin(m, data.compBlock
                                                     -> off0RhoBasisH2)
                       + kp(Id(m * d), hPrimes.second);           // superblock
    return lanczos(hSuper, psiGround, data.lancTolerance);
};

MatrixX_t TheBlock::createPrimeToRhoBasis(const MatrixX_t& rho, int mMax,
                                          double& cumulativeTruncationError)
{
    SelfAdjointEigenSolver<MatrixX_t> rhoSolver(rho);
                                             // find density matrix eigenstates
    int md = m * d,
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

MatrixX_t TheBlock::changeBasis(const MatrixX_t& mat, const TheBlock* block)
                                const
{
    return block -> primeToRhoBasis.adjoint() * mat * block -> primeToRhoBasis;
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            rmMatrixX_t& psiGround, int skips)
                                            const
{
    matPair hPrimes = createHprimes(data);
                                      // expanded system and environment blocks
    double gsEnergy = solveHSuper(hPrimes, data, psiGround);
                                                      // calculate ground state
    return FinalSuperblock(gsEnergy, data.ham.lSys, psiGround, m,
                           data.compBlock -> m, skips);
};

obsMatrixX_t TheBlock::obsChangeBasis(const obsMatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
