#include "FinalSuperblock.h"

#define obsId_d Matrix<obsScalar, d, d>::Identity()
#ifndef differentScalars
    #define obsProjectBlock projectBlock
    #define obsProjectBASCoupling projectBASCoupling
    #define obsProjectFreeSite projectFreeSite
#endif

using namespace Eigen;

FinalSuperblock::FinalSuperblock(double gsEnergy, int lSup,
                                 const rmMatrixX_t& psiGround, int mS, int mE,
                                 int skips)
    : gsEnergy(gsEnergy), lSup(lSup), psiGround(psiGround), mS(mS), mE(mE),
      skips(skips)
{
    if(lSup % 2)
    {
        lS = (lSup - 1)/2;
        lE = (lSup - 3)/2;
    }
    else
        lS = lE = lSup / 2 - 1;
};

double FinalSuperblock::expValue(const opsVec& ops,
                                 std::vector<TheBlock>& westBlocks,
                                 std::vector<TheBlock>& eastBlocks)
{
    sysBlockOps.clear();
    envBlockOps.clear();
    lFreeSite = obsId_d;
    rFreeSite = obsId_d;
    bool opInSysBlock = false,
         opAtlFreeSite = false,
         opAtrFreeSite = false,
         opInEnvBlock = false;
    for(const auto& opToEvaluate : ops) // split up operators into their blocks
    {
        int site = opToEvaluate.second;
        if(site < lS)
        {
            opInSysBlock = true;
            placeOp(opToEvaluate, sysBlockOps, true);
        }
        else if(site == lS)     // if necessary, assign left free-site operator
        {
            opAtlFreeSite = true;
            lFreeSite *= opToEvaluate.first;
        }
        else if(site == lS + 1)                 // and right free-site operator
        {
            opAtrFreeSite = true;
            rFreeSite *= opToEvaluate.first;
        }
        else
        {
            opInEnvBlock = true;
            placeOp(opToEvaluate, envBlockOps, false);
        };
    };
    psiGround.resize(mS * d * mE * d, 1);
    switch(8 * opInSysBlock + 4 * opAtlFreeSite + 2 * opInEnvBlock
           + opAtrFreeSite)                              // evaluate observable
    {
        case 0:                                                 // no operators
            return 1.;
        case 1:                                   // op at right-hand free site
            return obsRe(psiGround.col(0).eval().dot(actRFreeSite()));
        case 2:                                     // ops in environment block
            return obsRe(psiGround.col(0).eval().dot(actEnvBlock(eastBlocks)));
        case 3:            // ops in environment block and right-hand free site
            return obsRe(psiGround.col(0).eval()
                         .dot(actEnvBlockRFreeSite(eastBlocks)));
        case 4:                                         // op at left free site
            return obsRe(psiGround.col(0).eval().dot(actLFreeSite()));
        case 5:                                       // ops at both free sites
            return obsRe(actLFreeSite().dot(actRFreeSite()));
        case 6:          // ops in environment block and at left-hand free site
            return obsRe(actEnvBlock(eastBlocks).dot(actLFreeSite()));
        case 7:                 // ops at both free sites and environment block
            return obsRe(actLFreeSite().dot(actEnvBlockRFreeSite(eastBlocks)));
        case 8:                                          // ops in system block
            return obsRe(psiGround.col(0).eval().dot(actSysBlock(westBlocks)));
        case 9:              // ops in system block and at right-hand free site
            return obsRe(actRFreeSite().dot(actSysBlock(westBlocks)));
        case 10:                        // ops in system and environment blocks
            return obsRe(actSysBlock(westBlocks).dot(actEnvBlock(eastBlocks)));
        case 11: // ops in system and environment blocks and right-hand free site
            return obsRe(actSysBlock(westBlocks)
                         .dot(actEnvBlockRFreeSite(eastBlocks)));
        case 12:                 // ops in system block and left-hand free site
            return obsRe(psiGround.col(0).eval()
                         .dot(actSysBlockLFreeSite(westBlocks)));
        case 13:                     // ops in system block and both free sites
            return obsRe(actSysBlockLFreeSite(westBlocks)
                         .dot(actRFreeSite()));
        case 14: // ops in system and environment blocks and left-hand free site
            return obsRe(actSysBlockLFreeSite(westBlocks)
                         .dot(actEnvBlock(eastBlocks)));
        case 15:    // ops in system and environment blocks and both free sites
            return obsRe(actSysBlockLFreeSite(westBlocks)
                         .dot(actEnvBlockRFreeSite(eastBlocks)));
        default:
            return 0.;
        // can't actually reach this case - just here to quiet compiler warning
    };
};

void FinalSuperblock::placeOp(const std::pair<obsMatrixD_t, int>& op,
                              opsMap& blockSide, bool systemSide)
{
    int lhSite = (systemSide ? op.second : lSup - 1 - op.second);
    if(blockSide.count(lhSite))      // already an observable at this site?
        blockSide[lhSite] *= op.first;
    else
        blockSide.insert(std::pair<int, obsMatrixD_t>(lhSite, op.first));
};

obsVectorX_t FinalSuperblock::actSysBlock(std::vector<TheBlock>& westBlocks)
{
    psiGround.resize(mS, d * mE * d);
    obsMatrixX_t oPsi = rhoBasisRep(sysBlockOps, westBlocks, lS) * psiGround;
    oPsi.resize(mS * d * mE * d, 1);
    return oPsi;
};

obsVectorX_t FinalSuperblock::actLFreeSite()
{
    int mSd = mS * d,
        mEd = mE * d;
    psiGround.resize(mSd, mEd);
    rmMatrixX_t psiT = psiGround.transpose();
    psiT.resize(mEd * mS, d);
    obsMatrixX_t psiTOT = psiT * lFreeSite.transpose();
    psiTOT.resize(mEd, mSd);
    psiTOT.transposeInPlace();
    psiTOT.resize(mSd * mEd, 1);
    return psiTOT;
};

obsVectorX_t FinalSuperblock::actSysBlockLFreeSite(std::vector<TheBlock>&
                                                   westBlocks)
{
    obsMatrixX_t oLPsi = actLFreeSite();
    oLPsi.resize(mS, d * mE * d);
    obsMatrixX_t oSOLPsi = rhoBasisRep(sysBlockOps, westBlocks, lS) * oLPsi;
    oSOLPsi.resize(mS * d * mE * d, 1);
    return oSOLPsi;
};

obsVectorX_t FinalSuperblock::actEnvBlock(std::vector<TheBlock>& eastBlocks)
{
    int mSd = mS * d;
    psiGround.resize(mSd, mE * d);
    rmMatrixX_t psiT = psiGround.transpose();
    psiT.resize(mE, d * mSd);
    obsMatrixX_t oDagPsiT
        = rhoBasisRep(envBlockOps, eastBlocks, lE).adjoint() * psiT;
    oDagPsiT.resize(mE * d, mSd);
    oDagPsiT.transposeInPlace();
    oDagPsiT.resize(mSd * mSd, 1);
    return oDagPsiT;
};

obsVectorX_t FinalSuperblock::actRFreeSite()
{
    psiGround.resize(mS * d * mE, d);
    obsMatrixX_t psiOStar = psiGround * rFreeSite.conjugate();
    psiOStar.resize(mS * d * mE * d, 1);
    return psiOStar;
};

obsVectorX_t FinalSuperblock::actEnvBlockRFreeSite(std::vector<TheBlock>&
                                                   eastBlocks)
{
    obsMatrixX_t psiOEDag = actEnvBlock(eastBlocks);
    psiOEDag.resize(mS * d * mE, d);
    obsMatrixX_t psiOEDagORStar = psiOEDag * rFreeSite.conjugate();
    psiOEDagORStar.resize(mS * d * mE * d, 1);
    return psiOEDagORStar;
};

obsMatrixX_t FinalSuperblock::rhoBasisRep(const opsMap& blockOps,
                                          std::vector<TheBlock>& blocks,
                                          int blockSize) const
{
    obsMatrixX_t rhoBasisBlockOp;
    const auto firstOp = blockOps.begin(),    // first nontrivial site operator
               opsEnd = blockOps.end();                             // last one
    auto currentOp = firstOp;
    int currentSite;
    if(currentOp -> first < skips)                  // if necessary, do edge ED
    {
        if(currentOp -> first == 0)                 // at first site in system?
        {
            currentSite = 0;
            rhoBasisBlockOp = currentOp++ -> second;
        }
        else
        {
            currentSite = firstOp -> first - 1;
                              // move to last consecutive site with no operator
            rhoBasisBlockOp
                = obsMatrixX_t::Identity(blocks[currentSite].blockParts.m,
                                         blocks[currentSite].blockParts.m);
        };
        for(; currentSite < skips; currentSite++)
                         // perform ED tensoring until you reach end of edge ED
            rhoBasisBlockOp
                = kp(rhoBasisBlockOp, 
                     (   currentOp != opsEnd
                      && currentSite == currentOp -> first - 1 ?
                      currentOp++ -> second :
                      obsId_d)).eval();
    }
    else if(currentOp -> first == skips)
    {
        currentSite = skips;
        rhoBasisBlockOp
            = kp(obsMatrixX_t::Identity(blocks[currentSite - 1].blockParts.m,
                                        blocks[currentSite - 1].blockParts.m),
                 currentOp++ -> second);
    }
    else
    {
        currentSite = currentOp -> first - 1;
        rhoBasisBlockOp
            = blocks[currentSite++].obsProjectFreeSite(currentOp++ -> second);
    };
    for(int end = blockSize - 1; currentSite < end; currentSite++)
                                      // project rhoBasisBlockOp into new bases
        rhoBasisBlockOp
            = (currentOp != opsEnd && currentSite == currentOp -> first - 1 ?
               blocks[currentSite].obsProjectBASCoupling(rhoBasisBlockOp,
                                                         currentOp++ -> second) :
               blocks[currentSite].obsProjectBlock(rhoBasisBlockOp));
    return rhoBasisBlockOp;
};
