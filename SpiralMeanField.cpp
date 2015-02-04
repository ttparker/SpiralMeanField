#include "TheBlock.h"

#define j1 couplingConstants[0]
#define j2 couplingConstants[1]
#define jprime couplingConstants[2]
#define hIn couplingConstantsIn[3]
#define sigmaplus siteBasisH2[0]
#define sigmaz siteBasisH2[1]
#define sigmaminus siteBasisH2[2]
#define rhoBasisSigmaplus rhoBasisH2[0]
#define rhoBasisSigmaz rhoBasisH2[1]
#define off0RhoBasisSigmaplus off0RhoBasisH2[0]
#define off0RhoBasisSigmaz off0RhoBasisH2[1]
#define off1RhoBasisSigmaplus off1RhoBasisH2[0]
#define off1RhoBasisSigmaz off1RhoBasisH2[1]

#define lS blockParts.l
#define mS blockParts.m
#define lE compBlockParts.l
#define mE compBlockParts.m

using namespace Eigen;

Hamiltonian::Hamiltonian()
{
    siteBasisH2.resize(nCouplingOperators);
    sigmaplus << 0., 1.,
                 0., 0.;                           // define one-site operators
    sigmaminus << 0., 0.,
                  1., 0.;
    sigmaz << 1.,  0.,
              0., -1.;
    sigmax << 0., 1.,
              1., 0.;
    sigmay <<                          0., std::complex<double>(0., -1.),
              std::complex<double>(0., 1),                            0.;
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
    MatrixX_t plusMinus = kp(rhoBasisSigmaplus, sigmaminus);
    return couplingConstants[j - 1] * (  kp(rhoBasisSigmaz, sigmaz)
                                       + 2 * (plusMinus + plusMinus.adjoint()));
};

MatrixD_t Hamiltonian::h1(int sitesFromWestEnd) const
{
    return -(   hFromIntSpins(sitesFromWestEnd, 0)      * sigmax
             +  hFromIntSpins(sitesFromWestEnd, 1)      * sigmay
             + (hFromIntSpins(sitesFromWestEnd, 2) + h) * sigmaz);
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
        = j2 * (  sysBlockRSiteCoupling(blockParts.off0RhoBasisSigmaz, sigmaz,
                                        x, mS, mE)
                + 2 * (  sysBlockRSiteCoupling(blockParts.off0RhoBasisSigmaplus,
                                               sigmaminus, x, mS, mE)
                       + sysBlockRSiteCoupling(blockParts.off0RhoBasisSigmaplus
                                               .adjoint(),
                                               sigmaplus, x, mS, mE)));
    y1234.resize(totalDimension, 1);
    
    // act with left-hand side of system:
    rmMatrixX_t y2341 = (blockParts.hS * x).transpose(), // act with system block
                cycledx = x.transpose();
    y2341.resize(d, mEd * mS);
    y2341.noalias()
        +=  j1 * (  sysBlockLSiteCoupling(blockParts.off0RhoBasisSigmaz, sigmaz,
                                          cycledx, mS, mE)
                  + 2 * (  sysBlockLSiteCoupling(blockParts.off0RhoBasisSigmaplus,
                                                 sigmaminus, cycledx, mS, mE)
                         + sysBlockLSiteCoupling(blockParts.off0RhoBasisSigmaplus
                                                 .adjoint(),
                                                 sigmaplus, cycledx, mS, mE)))
          + j2 * (  sysBlockLSiteCoupling(blockParts.off1RhoBasisSigmaz, sigmaz,
                                          cycledx, mS, mE)
                  + 2 * (  sysBlockLSiteCoupling(blockParts.off1RhoBasisSigmaplus,
                                                 sigmaminus, cycledx, mS, mE)
                         + sysBlockLSiteCoupling(blockParts.off1RhoBasisSigmaplus
                                                 .adjoint(),
                                                 sigmaplus, cycledx, mS, mE)));
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
            = j1 * (  freeSiteCouplingTerm(sigmaz, sigmaz, RHS, mE)
                    + 2 * (  freeSiteCouplingTerm(sigmaplus, sigmaminus, RHS, mE)
                           + freeSiteCouplingTerm(sigmaminus, sigmaplus, RHS, mE)));
    };
    yFSBond.resize(totalDimension, 1);
    
    // act with environment-block-left-hand-free-site bond:
    x.resize(mSd, mEd);
    cycledx = x.transpose();
    cycledx.resize(mE, d * mSd);
    rmMatrixX_t y3412
        = j2 * (  envBlockLSiteCoupling(compBlockParts.off0RhoBasisSigmaz, sigmaz,
                                        cycledx, mS, mE)
                + 2 * (  envBlockLSiteCoupling(compBlockParts
                                               .off0RhoBasisSigmaplus,
                                               sigmaminus, cycledx, mS, mE)
                       + envBlockLSiteCoupling(compBlockParts
                                               .off0RhoBasisSigmaplus.adjoint(),
                                               sigmaplus, cycledx, mS, mE)));
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
        +=  j1 * (  envBlockRSiteCoupling(compBlockParts.off0RhoBasisSigmaz,
                                          sigmaz, cycledx, mS, mE)
                  + 2 * (  envBlockRSiteCoupling(compBlockParts
                                                 .off0RhoBasisSigmaplus,
                                                 sigmaminus, cycledx, mS, mE)
                         + envBlockRSiteCoupling(compBlockParts
                                                 .off0RhoBasisSigmaplus.adjoint(),
                                                 sigmaplus, cycledx, mS, mE)))
          + j2 * (  envBlockRSiteCoupling(compBlockParts.off1RhoBasisSigmaz,
                                          sigmaz, cycledx, mS, mE)
                  + 2 * (  envBlockRSiteCoupling(compBlockParts
                                                 .off1RhoBasisSigmaplus,
                                                 sigmaminus, cycledx, mS, mE)
                         + envBlockRSiteCoupling(compBlockParts
                                                 .off1RhoBasisSigmaplus.adjoint(),
                                                 sigmaplus, cycledx, mS, mE)));
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
                    -> projectBASCoupling(block -> blockParts.off0RhoBasisSigmaz,
                                          sigmaz)
                  + 2 * (  block
                           -> projectBASCoupling(block
                                                 -> blockParts
                                                 .off0RhoBasisSigmaplus,
                                                 sigmaminus)
                         + block
                           -> projectBASCoupling(block
                                                 -> blockParts
                                                 .off0RhoBasisSigmaplus.adjoint(),
                                                 sigmaplus)))
          + j2 * (  block
                    -> projectBASCoupling(block -> blockParts.off1RhoBasisSigmaz,
                                          sigmaz)
                  + 2 * (  block
                           -> projectBASCoupling(block
                                                 -> blockParts
                                                 .off1RhoBasisSigmaplus,
                                                 sigmaminus)
                         + block
                           -> projectBASCoupling(block
                                                 -> blockParts
                                                 .off1RhoBasisSigmaplus.adjoint(),
                                                 sigmaplus)));
};
