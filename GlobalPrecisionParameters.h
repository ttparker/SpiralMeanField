#ifndef GPP_H
#define GPP_H

const int globalMinLancIters = 4,
          globalMaxLancIters = 100;
const double fallbackLancTolerance = 1.e-4,
             // reduced error tolerance to accept if Lanczos fails to converge
             degenerateDMCutoff = 1.e-4,
             // smallest relative difference between consecutive density-matrix
             // eigenvalues between which to truncate the Hamiltonian projection
             observableThreshold = 1.e-8;
             // if any observable's absolute value is smaller than this, then
             // round it to zero

#endif
