#ifndef FINALSUPERBLOCK_H
#define FINALSUPERBLOCK_H

#include <map>
#include "TheBlock.h"

typedef std::vector<std::pair<obsMatrixD_t, int>,
                    Eigen::aligned_allocator<std::pair<obsMatrixD_t, int>>>
                    opsVec;
typedef std::map<int, obsMatrixD_t, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, obsMatrixD_t>>>
                 opsMap;

class FinalSuperblock
{
    public:
        double gsEnergy;                         // returns ground-state energy
        
        FinalSuperblock(double gsEnergy, int lSupFinal,
                        const rmMatrixX_t& psiGround,
                        int mSFinal, int mEFinal, int skips);
        double expValue(const opsVec& ops, std::vector<TheBlock>& westBlocks,
                        std::vector<TheBlock>& eastBlocks);
        // calculates exectation value of a combination of single-site operators
    
    private:
        int lSupFinal;                                     // final system size
        rmMatrixX_t psiGround;                 // final superblock ground state
        int lSFinal,                        // final length of the system block
            lEFinal,                   // final length of the environment block
            mSFinal,           // final number of states stored in system block
            mEFinal,      // final number of states stored in environment block
            skips;
        
        void placeOp(const std::pair<obsMatrixD_t, int>& op, opsMap& blockSide,
                     bool systemSide);
                    // assign each one-site observable to the appropriate block
        obsMatrixX_t rhoBasisRep(const opsMap& blockOps,
                                 std::vector<TheBlock>& blocks, int blockSize)
                                 const;
                // converts single-site operators into the system block basis
};

#endif
