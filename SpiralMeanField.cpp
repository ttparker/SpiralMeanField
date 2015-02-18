#include "TheBlock.h"

#define j1 couplingConstants[0]
#define j2 couplingConstants[1]
#define jprime couplingConstants[2]
#define hIn couplingConstantsIn[3]
#define splus siteBasisH2[0]
#define sz siteBasisH2[1]
#define sminus siteBasisH2[2]
#define rhoBasissplus rhoBasisH2[0]
#define rhoBasissz rhoBasisH2[1]
#define off0RhoBasissplus off0RhoBasisH2[0]
#define off0RhoBasissz off0RhoBasisH2[1]
#define off1RhoBasissplus off1RhoBasisH2[0]
#define off1RhoBasissz off1RhoBasisH2[1]

#define lS blockParts.l
#define mS blockParts.m
#define lE compBlockParts.l
#define mE compBlockParts.m

using namespace Eigen;

Hamiltonian::Hamiltonian()
{
    siteBasisH2.resize(nCouplingOperators);
    splus << 0., 1.,
             0., 0.;                               // define one-site operators
    sminus << 0., 0.,
              1., 0.;
    sz << .5,  0.,
          0., -.5;
    sx << 0., .5,
          .5, 0.;
    sy <<                           0., std::complex<double>(0., -.5),
          std::complex<double>(0., .5),                            0.;
};

void Hamiltonian::setParams(const std::vector<double>& couplingConstantsIn,
                            int lSysIn)
{
    couplingConstants.assign(couplingConstantsIn.begin(),
                             couplingConstantsIn.end() - 1);
    lSys = lSysIn;
    h = hIn;
};

void Hamiltonian::calcEffectiveH(const MatrixXd& intSpins)
{
    hFromIntSpins = MatrixXd(lSys, 3);
    hFromIntSpins.row(0) = -jprime * intSpins.row(0);
                       // interstitial spins have hard-wall boundary conditions
    hFromIntSpins.bottomRows(lSys - 1)
        = -jprime * (intSpins.bottomRows(lSys - 1) + intSpins.topRows(lSys - 1));
};

rmMatrixX_t
    Hamiltonian::EDBASCoupling(int j, const std::vector<MatrixX_t>& rhoBasisH2)
    const
{
    MatrixX_t plusMinus = kp(rhoBasissplus, sminus);
    return couplingConstants[j - 1] * (  kp(rhoBasissz, sz)
                                       + (plusMinus + plusMinus.adjoint()) / 2);
};

MatrixD_t Hamiltonian::h1(int sitesFromWestEnd) const
{
    return -(   hFromIntSpins(sitesFromWestEnd, 0)      * sx
             +  hFromIntSpins(sitesFromWestEnd, 1)      * sy
             + (hFromIntSpins(sitesFromWestEnd, 2) + h) * sz);
};

Matrix<hamScalar, Dynamic, 1>
    Hamiltonian::act(const effectiveHams& blockParts,
                     const effectiveHams& compBlockParts, rmMatrixX_t x,
                     bool sweepingEast) const
{
    int mSd = mS * d,
        mEd = mE * d,
        totalDimension = mSd * mEd,
                                // dimension of state on which Hamiltonian acts
        lFreeSiteDistFromWestEnd,
        rFreeSiteDistFromWestEnd;
    if(sweepingEast)
    {
        lFreeSiteDistFromWestEnd = lS + 1;
        rFreeSiteDistFromWestEnd = lSys - 2 - lE;
    }
    else
    {
        lFreeSiteDistFromWestEnd = lSys - 2 - lS;
        rFreeSiteDistFromWestEnd = lE + 1;
    };
    
    // act with system-block-right-hand-free-site bond:
    x.resize(mS, d * mEd);
    rmMatrixX_t y1234
        = j2 * (  sysBlockRSiteCoupling(blockParts.off0RhoBasissz, sz, x, mS, mE)
                + (  sysBlockRSiteCoupling(blockParts.off0RhoBasissplus, sminus,
                                           x, mS, mE)
                   + sysBlockRSiteCoupling(blockParts.off0RhoBasissplus
                                           .adjoint(), splus, x, mS, mE)) / 2);
    y1234.resize(totalDimension, 1);
    
    // act with left-hand side of system:
    rmMatrixX_t y2341 = (blockParts.hS * x).transpose(), // act with system block
                cycledx = x.transpose();
    y2341.resize(d, mEd * mS);
    y2341.noalias()
        +=  j1 * (  sysBlockLSiteCoupling(blockParts.off0RhoBasissz, sz,
                                          cycledx, mS, mE)
                  + (  sysBlockLSiteCoupling(blockParts.off0RhoBasissplus,
                                             sminus, cycledx, mS, mE)
                     + sysBlockLSiteCoupling(blockParts.off0RhoBasissplus
                                             .adjoint(), splus, cycledx, mS, mE))
                    / 2)
          + j2 * (  sysBlockLSiteCoupling(blockParts.off1RhoBasissz, sz,
                                          cycledx, mS, mE)
                  + (  sysBlockLSiteCoupling(blockParts.off1RhoBasissplus,
                                             sminus, cycledx, mS, mE)
                     + sysBlockLSiteCoupling(blockParts.off1RhoBasissplus
                                             .adjoint(),
                                             splus, cycledx, mS, mE)) / 2);
                             // act with system-block-left-hand-free-site bonds
    cycledx.resize(d, mEd * mS);
    y2341.noalias() += h1(lFreeSiteDistFromWestEnd) * cycledx;
                                                // act with left-hand free site
    y2341.resize(d * mEd, mS);
    y2341.transposeInPlace();
    y2341.resize(totalDimension, 1);
    
    // act with free sites' bond:
    rmMatrixX_t yFSBond(mS, d * mEd);
    for(int i = 0; i < mS; i++)
    {
        rmMatrixX_t RHS = x.row(i);
        RHS.resize(d, mEd);
        yFSBond.row(i)
            = j1 * (  freeSiteCouplingTerm(sz, sz, RHS, mE)
                    + (  freeSiteCouplingTerm(splus, sminus, RHS, mE)
                       + freeSiteCouplingTerm(sminus, splus, RHS, mE)) / 2);
    };
    yFSBond.resize(totalDimension, 1);
    
    // act with environment-block-left-hand-free-site bond:
    x.resize(mSd, mEd);
    cycledx = x.transpose();
    cycledx.resize(mE, d * mSd);
    rmMatrixX_t y3412
        = j2 * (  envBlockLSiteCoupling(compBlockParts.off0RhoBasissz, sz,
                                        cycledx, mS, mE)
                + (  envBlockLSiteCoupling(compBlockParts.off0RhoBasissplus,
                                           sminus, cycledx, mS, mE)
                   + envBlockLSiteCoupling(compBlockParts.off0RhoBasissplus
                                           .adjoint(),
                                           splus, cycledx, mS, mE)) / 2);
    y3412.resize(mEd, mSd);
    y3412.transposeInPlace();
    y3412.resize(totalDimension, 1);
    
    // act with right-hand half of system:
    x.resize(mSd * mE, d);
    cycledx = x.transpose();
    rmMatrixX_t y4123 = h1(rFreeSiteDistFromWestEnd) * cycledx;
                                               // act with right-hand free site
    y4123.resize(d * mSd, mE);
    y4123.noalias()
        +=  j1 * (  envBlockRSiteCoupling(compBlockParts.off0RhoBasissz, sz,
                                          cycledx, mS, mE)
                  + (  envBlockRSiteCoupling(compBlockParts.off0RhoBasissplus,
                                             sminus, cycledx, mS, mE)
                     + envBlockRSiteCoupling(compBlockParts.off0RhoBasissplus
                                             .adjoint(),
                                             splus, cycledx, mS, mE)) / 2)
          + j2 * (  envBlockRSiteCoupling(compBlockParts.off1RhoBasissz, sz,
                                          cycledx, mS, mE)
                  + (  envBlockRSiteCoupling(compBlockParts.off1RhoBasissplus,
                                             sminus, cycledx, mS, mE)
                     + envBlockRSiteCoupling(compBlockParts.off1RhoBasissplus
                                             .adjoint(),
                                             splus, cycledx, mS, mE)) / 2);
                       // act with environment-block-right-hand-free-site bonds
    cycledx.resize(d * mSd, mE);
    y4123.noalias() += cycledx * compBlockParts.hS.transpose();
                                                  // act with environment block
    y4123.resize(d, mSd * mE);
    y4123.transposeInPlace();
    y4123.resize(totalDimension, 1);
    
    return y1234 + y2341 + yFSBond + y3412 + y4123;
};

rmMatrixX_t Hamiltonian::sysBlockRSiteCoupling(const MatrixX_t& hSBond,
                                               const MatrixD_t& hSiteBond,
                                               const rmMatrixX_t& x, int mSIn,
                                               int mEIn) const
{
    rmMatrixX_t hBondx = hSBond * x;
    hBondx.resize(mSIn * d * mEIn, d);
    return hBondx * hSiteBond.transpose();
};

rmMatrixX_t Hamiltonian::sysBlockLSiteCoupling(const MatrixX_t& hSBond,
                                               const MatrixD_t& hSiteBond,
                                               const rmMatrixX_t& cycledx,
                                               int mSIn, int mEIn) const
{
    rmMatrixX_t hBondx = cycledx * hSBond.transpose();
                                     // act with system-block half of bond term
    hBondx.resize(d, mEIn * d * mSIn);
    return hSiteBond * hBondx; // act with left-hand-free-site half of bond term
};

Matrix<hamScalar, 1, Dynamic>
    Hamiltonian::freeSiteCouplingTerm(const MatrixD_t& lFreeSiteBond,
                                      const MatrixD_t& rFreeSiteBond,
                                      const rmMatrixX_t& RHS, int mEIn) const
{
    rmMatrixX_t hSiteBondRHS = lFreeSiteBond * RHS;
                                   // act with left-hand-free-site half of bond
    hSiteBondRHS.resize(mEIn * d, d);
    rmMatrixX_t hSiteBondRHShSiteBond
        = hSiteBondRHS * rFreeSiteBond.transpose();
                                  // act with right-hand-free-site half of bond
    hSiteBondRHShSiteBond.resize(1, d * mEIn * d);
    return hSiteBondRHShSiteBond;
};

rmMatrixX_t Hamiltonian::envBlockLSiteCoupling(const MatrixX_t& hEBond,
                                               const MatrixD_t& hSiteBond,
                                               const rmMatrixX_t& cycledx,
                                               int mSIn, int mEIn) const
{
    rmMatrixX_t hBondx = hEBond * cycledx;
    hBondx.resize(mEIn * d * mSIn, d);
    return hBondx * hSiteBond.transpose();
};

rmMatrixX_t Hamiltonian::envBlockRSiteCoupling(const MatrixX_t& hEBond,
                                               const MatrixD_t& hSiteBond,
                                               const rmMatrixX_t& cycledx,
                                               int mSIn, int mEIn) const
{
    rmMatrixX_t hBondx = hSiteBond * cycledx;
                             // act with right-hand free site half of bond term
    hBondx.resize(d * mSIn * d, mEIn);
    return hBondx * hEBond.transpose();
                                // act with environment-block half of bond term
};

rmMatrixX_t Hamiltonian::projectedBASCouplings(TheBlock * block) const
{
    return  j1 * (  block
                    -> projectBASCoupling(block -> blockParts.off0RhoBasissz, sz)
                  + (  block
                       -> projectBASCoupling(block -> blockParts.off0RhoBasissplus,
                                             sminus)
                     + block
                       -> projectBASCoupling(block -> blockParts
                                             .off0RhoBasissplus.adjoint(),
                                             splus)) / 2)
          + j2 * (  block
                    -> projectBASCoupling(block -> blockParts.off1RhoBasissz, sz)
                  + (  block
                       -> projectBASCoupling(block -> blockParts.off1RhoBasissplus,
                                             sminus)
                     + block
                       -> projectBASCoupling(block -> blockParts.off1RhoBasissplus
                                                      .adjoint(),
                                             splus)) / 2);
};
