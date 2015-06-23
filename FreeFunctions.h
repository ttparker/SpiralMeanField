#ifndef FREEFUNCTIONS_H
#define FREEFUNCTIONS_H

#include "FinalSuperblock.h"

rmMatrixX_t randomSeed(const TheBlock& leftBlock, const TheBlock& rightBlock);
                                            // outputs random normalized vector
void reflectPredictedPsi(rmMatrixX_t& psiGround, const TheBlock& bigBlock,
                         const TheBlock& littleBlock);
                                               // when you reach edge of system
Eigen::VectorXd oneSiteExpValues(const obsMatrixD_t& oneSiteOp,
                                 int rangeOfObservables, int lSys,
                                 FinalSuperblock& hSuperFinal,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks);
Eigen::MatrixXd twoSiteExpValues(const obsMatrixD_t& firstTwoSiteOp,
                                 const obsMatrixD_t& secondTwoSiteOp,
                                 int rangeOfObservables, int lSys,
                                 FinalSuperblock& hSuperFinal,
                                 std::vector<TheBlock>& westBlocks,
                                 std::vector<TheBlock>& eastBlocks,
                                 std::ofstream& fileout);

#endif
