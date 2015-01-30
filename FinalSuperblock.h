#ifndef FINALSUPERBLOCK_H
#define FINALSUPERBLOCK_H

#include <map>
#include "TheBlock.h"

typedef std::map<int, obsMatrixD_t, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, obsMatrixD_t>>>
                 opsMap;
typedef Eigen::Matrix<obsScalar, Eigen::Dynamic, 1> obsVectorX_t;

class FinalSuperblock
{
    public:
        double gsEnergy;                         // returns ground-state energy
        
        FinalSuperblock(double gsEnergy, int lSupFinal,
                        const rmMatrixX_t& psiGround,
                        int mSFinal, int mEFinal, int skips);
        double expValue(const opsVec& ops, std::vector<TheBlock>& westBlocks,
                        std::vector<TheBlock>& eastBlocks);
      // calculates expectation value of a combination of single-site operators
    
    private:
        // system information set upon initialization:
        int lSup;                                          // final system size
        rmMatrixX_t psiGround;                 // final superblock ground state
        int lS,                             // final length of the system block
            lE,                        // final length of the environment block
            mS,                // final number of states stored in system block
            mE,           // final number of states stored in environment block
            skips;
        // observables that are currently being evaluated:
        opsMap sysBlockOps,
                      // observable operators that will act on the system block
               envBlockOps;                       // same for environment block
        obsMatrixD_t lFreeSite,
                     rFreeSite;

        void placeOp(const std::pair<obsMatrixD_t, int>& op, opsMap& blockSide,
                     bool systemSide);
                    // assign each one-site observable to the appropriate block
        obsVectorX_t actSysBlock(std::vector<TheBlock>& leftBlocks),
                                    // act the system block on the ground state
                     actLFreeSite(),
                            // acts the left-hand free site on the ground state
                     actSysBlockLFreeSite(std::vector<TheBlock>& leftBlocks),
                                     // acts both the system block and the
                                     // left-hand free site on the ground state
                     actEnvBlock(std::vector<TheBlock>& rightBlocks),
               // acts the ADJOINT of the environment block on the ground state
                     actRFreeSite(),      // acts the ADJOINT of the right-hand
                                          // free site on the ground state
                     actEnvBlockRFreeSite(std::vector<TheBlock>& rightBlocks);
                    // acts both the ADJOINT of the right-hand free site and
                    // the ADJOINT of the environment block on the ground state
        obsMatrixX_t rhoBasisRep(const opsMap& blockOps,
                                 std::vector<TheBlock>& blocks, int blockSize)
                                 const;
                // converts single-site operators into the system block basis
};

#endif
