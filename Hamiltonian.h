#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"

#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)
#define Id_d Matrix<double, d, d>::Identity()       // one-site identity matrix

typedef std::vector<MatrixD_t, Eigen::aligned_allocator<MatrixD_t>> vecMatD_t;

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
        vecMatD_t siteBasisH2;       // site-basis coupling operators - the
                                     // independent ones should be listed first
        double h;                           // external magnetic field strength
        std::vector<double> cosList,       // cos(k i) for each i in the system
                            sinList;       // same for sin(k i)
        MatrixD_t sigmax;                             // for the external field
        
        MatrixX_t blockAdjacentSiteJoin(int j,
                                        const std::vector<MatrixX_t>& rhoBasisH2)
                                        const,
                                         // j gives the j-1th coupling constant
                  lBlockrSiteJoin(const std::vector<MatrixX_t>& off0RhoBasisH2,
                                  int compm) const,
                  lSiterBlockJoin(int m,
                                  const std::vector<MatrixX_t>& off0RhoBasisH2)
                                  const,
                  siteSiteJoin(int m, int compm) const;
                                           // joins the two free sites together
        MatrixD_t h1(int sitesFromWestEnd) const;
                        // returns the external magnetic field as a function of
                        // position, measured from the initially leftmost site
    
    friend class TheBlock;
};

#endif
