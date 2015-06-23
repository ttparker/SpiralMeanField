#include <fstream>
#include "FinalSuperblock.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

rmMatrixX_t randomSeed(const TheBlock& leftBlock, const TheBlock& rightBlock)
{
    return rmMatrixX_t::Random(leftBlock.blockParts.m * d
                               * rightBlock.blockParts.m * d, 1).normalized();
};

void reflectPredictedPsi(rmMatrixX_t& psiGround, const TheBlock& bigBlock,
                         const TheBlock& littleBlock)
{
    psiGround.resize(bigBlock.blockParts.m * d, littleBlock.blockParts.m * d);
    psiGround.transposeInPlace();
    psiGround.resize(bigBlock.blockParts.m * d * littleBlock.blockParts.m * d, 1);
};

VectorXd oneSiteExpValues(const obsMatrixD_t& oneSiteOp, int lSys,
                          FinalSuperblock& hSuperFinal,
                          std::vector<TheBlock>& westBlocks,
                          std::vector<TheBlock>& eastBlocks)
{
    opsVec ops;                     // list of observable single-site operators
    ops.push_back(std::make_pair(oneSiteOp, 0));
    VectorXd oneSiteVals(lSys);
    for(int i = 0; i < lSys; i++)
    {
        ops[0].second = i;
        double exactValue = hSuperFinal.expValue(ops, westBlocks, eastBlocks);
        oneSiteVals(i) = std::abs(exactValue) < observableThreshold ?
                         0. : exactValue;
    };
    return oneSiteVals;
};

double expiDotj(int i, int j, FinalSuperblock& hSuperFinal,
                std::vector<TheBlock>& leftBlocks,
                std::vector<TheBlock>& rightBlocks)
{
    obsMatrixD_t splus,
                 sminus,
                 sz;
    splus << 0, 1,
             0, 0;
    sminus << 0, 0,
              1, 0;
    sz << .5,   0,
           0, -.5;
    opsVec ops;                     // list of observable single-site operators
    ops.reserve(2);
    ops.push_back(std::make_pair(splus, i));
    ops.push_back(std::make_pair(sminus, j));
    double dotExpValue = hSuperFinal.expValue(ops, leftBlocks, rightBlocks) / 2;
    ops[0].first = sminus;
    ops[1].first = splus;
    dotExpValue += hSuperFinal.expValue(ops, leftBlocks, rightBlocks) / 2;
    ops[0].first = sz;
    ops[1].first = sz;
    dotExpValue += hSuperFinal.expValue(ops, leftBlocks, rightBlocks);
    return dotExpValue;
};
