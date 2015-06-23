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

VectorXd oneSiteExpValues(const obsMatrixD_t& oneSiteOp, int rangeOfObservables,
                          int lSys, FinalSuperblock& hSuperFinal,
                          std::vector<TheBlock>& leftBlocks,
                          std::vector<TheBlock>& rightBlocks)
{
    opsVec ops;                     // list of observable single-site operators
    ops.push_back(std::make_pair(oneSiteOp, 0));
    VectorXd oneSiteVals(rangeOfObservables);
    int start = (lSys - rangeOfObservables) / 2;
    for(int i = 0; i < rangeOfObservables; i++)
    {
        ops[0].second = start + i;
        double exactValue = hSuperFinal.expValue(ops, leftBlocks, rightBlocks);
        oneSiteVals(i) = std::abs(exactValue) < observableThreshold ?
                         0. : exactValue;
    };
    return oneSiteVals;
};

MatrixXd twoSiteExpValues(const obsMatrixD_t& firstTwoSiteOp,
                          const obsMatrixD_t& secondTwoSiteOp,
                          int rangeOfObservables,
                          int lSys, FinalSuperblock& hSuperFinal,
                          std::vector<TheBlock>& westBlocks,
                          std::vector<TheBlock>& eastBlocks,
                          std::ofstream& fileout)
{
    opsVec ops;                     // list of observable single-site operators
    ops.reserve(2);
    ops.push_back(std::make_pair(firstTwoSiteOp, 0));
    ops.push_back(std::make_pair(secondTwoSiteOp, 0));
    MatrixXd correlationFunction = MatrixXd::Zero(rangeOfObservables,
                                                  rangeOfObservables);
    int start = (lSys - rangeOfObservables) / 2;
    for(int i = 0; i < rangeOfObservables; i++)
        for(int j = 0; j < rangeOfObservables; j++)
        {
            ops[0].second = start + i;
            ops[1].second = start + j;
            double exactValue = hSuperFinal.expValue(ops, westBlocks,
                                                     eastBlocks);
            correlationFunction(i, j) = std::abs(exactValue)
                                            < observableThreshold ?
                                        0. : exactValue;
        };
    fileout << "Two-site correlation function:\n" << correlationFunction
            << std::endl << std::endl;
    return correlationFunction;
};
