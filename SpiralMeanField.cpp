#include "Hamiltonian.h"

#define j1 couplingConstants[0]
#define j2 couplingConstants[1]
#define sigmaplus siteBasisH2[0]
#define sigmaz siteBasisH2[1]
#define sigmaminus siteBasisH2[2]
#define rhoBasisSigmaplus rhoBasisH2[0]
#define rhoBasisSigmaz rhoBasisH2[1]
#define off0RhoBasisSigmaplus off0RhoBasisH2[0]
#define off0RhoBasisSigmaz off0RhoBasisH2[1]
#define off1RhoBasisSigmaplus off1RhoBasisH2[0]
#define off1RhoBasisSigmaz off1RhoBasisH2[1]

using namespace Eigen;

Hamiltonian::Hamiltonian()
{
    siteBasisH2.resize(nCouplingOperators);
    sigmaplus << 0., 1.,
                 0., 0.;                  // define one-site operators matrices
    sigmaminus << 0., 0.,
                  1., 0.;
    sigmaz << 1.,  0.,
              0., -1.;
    sigmax << 0., 1.,
              1., 0.;
};

void Hamiltonian::setParams(const std::vector<double>& couplingConstantsIn,
                            int lSysIn, double k)
{
    couplingConstants.assign(couplingConstantsIn.begin(),
                             couplingConstantsIn.end() - 1);
    lSys = lSysIn;
    h = -2 * couplingConstantsIn[2] * cos(k / 2);
    cosList.clear();
    sinList.clear();
    cosList.reserve(lSysIn);
    sinList.reserve(lSysIn);
    for(int i = 0; i < lSysIn; i++)
    {
        cosList.push_back(cos(k * i));
        sinList.push_back(sin(k * i));
    };
}

MatrixX_t Hamiltonian::blockAdjacentSiteJoin(int j,
                                             const std::vector<MatrixX_t>& rhoBasisH2)
                                             const
{
    MatrixX_t plusMinus = kp(rhoBasisSigmaplus, sigmaminus);
    return couplingConstants[j - 1] *
        (kp(rhoBasisSigmaz, sigmaz) + 2 * (plusMinus + plusMinus.adjoint()));
};

MatrixX_t Hamiltonian::lBlockrSiteJoin(const std::vector<MatrixX_t>&
                                       off0RhoBasisH2, int compm) const
{
    MatrixX_t plusMinus = kp(kp(off0RhoBasisSigmaplus, Id(d * compm)),
                             sigmaminus);
    return j2 * (kp(kp(off0RhoBasisSigmaz, Id(d * compm)), sigmaz)
                 + 2 * (plusMinus + plusMinus.adjoint()));
};

MatrixX_t Hamiltonian::lSiterBlockJoin(int m,
                                       const std::vector<MatrixX_t>&
                                       off0RhoBasisH2) const
{
    MatrixX_t plusMinus = kp(sigmaplus, off0RhoBasisSigmaplus.adjoint());
    return j2 *
        kp(kp(Id(m), kp(sigmaz, off0RhoBasisSigmaz)
                     + 2 * (plusMinus + plusMinus.adjoint())),
           Id_d);
};

MatrixX_t Hamiltonian::siteSiteJoin(int m, int compm) const
{
    MatrixX_t plusMinus = kp(kp(sigmaplus, Id(compm)), sigmaminus);
    return j1 * kp(Id(m), kp(kp(sigmaz, Id(compm)), sigmaz)
                           + 2 * (plusMinus + plusMinus.adjoint()));
};

MatrixD_t Hamiltonian::h1(int sitesFromWestEnd) const
{
    return -h * (cosList[sitesFromWestEnd] * sigmax
                 + sinList[sitesFromWestEnd] * sigmaz);
};
