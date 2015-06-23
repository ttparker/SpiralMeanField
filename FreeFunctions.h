#ifndef FREEFUNCTIONS_H
#define FREEFUNCTIONS_H

#include "FinalSuperblock.h"

rmMatrixX_t randomSeed(const TheBlock& leftBlock, const TheBlock& rightBlock);
                                            // outputs random normalized vector
void reflectPredictedPsi(rmMatrixX_t& psiGround, const TheBlock& bigBlock,
                         const TheBlock& littleBlock);
                                               // when you reach edge of system
Eigen::VectorXd oneSiteExpValues(const obsMatrixD_t& oneSiteOp,
                                 int lSys, FinalSuperblock& hSuperFinal,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks);
double expiDotj(int i, int j, FinalSuperblock& hSuperFinal,
                std::vector<TheBlock>& leftBlocks,
                std::vector<TheBlock>& rightBlocks);
#endif
