#include <time.h>
#include <fstream>
#include "FreeFunctions.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

int main()
{
    clock_t start = clock();
    
    std::ifstream filein("Input/Input");
    if(!filein)
    {
        std::cerr << "Couldn't open input file." << std::endl;
        exit(EXIT_FAILURE);
    };
    
    // read in parameters that are constant across all trials:
    int nTrials;
    bool calcObservables;
    filein >> nTrials >> calcObservables;
    bool calcOneSiteExpValues;
    int indexOfOneSiteOp;
    obsMatrixD_t oneSiteOp;
    bool calcTwoSiteExpValues;
    int indexOfFirstTwoSiteOp,
        indexOfSecondTwoSiteOp;
    obsMatrixD_t firstTwoSiteOp,
                 secondTwoSiteOp;
    #include "ObservableOps.h"
    if(calcObservables)
    {
        filein >> calcOneSiteExpValues;
        if(calcOneSiteExpValues)
        {
            filein >> indexOfOneSiteOp;
            oneSiteOp = obsList[indexOfOneSiteOp];
        };
        filein >> calcTwoSiteExpValues;
        if(calcTwoSiteExpValues)
        {
            filein >> indexOfFirstTwoSiteOp >> indexOfSecondTwoSiteOp;
            firstTwoSiteOp = obsList[indexOfFirstTwoSiteOp];
            secondTwoSiteOp = obsList[indexOfSecondTwoSiteOp];
        };
    };
    
    std::ofstream infoout("Output/Info");
    if(infoout)
    {
        infoout << "Number of trials: " << nTrials
                << "\nCalculate observables? "
                << (calcObservables ? "Yes" : "No") << std::endl;
        if(calcObservables)
        {
            infoout << "Calculate one-site observables? ";
            if(calcOneSiteExpValues)
                infoout << "Yes\nIndex of one-site observable: "
                        << indexOfOneSiteOp << std::endl;
            else
                infoout << "No" << std::endl;
            infoout << "Calculate two-site observables? ";
            if(calcTwoSiteExpValues)
                infoout << "Yes\nIndices of two-site observables: "
                        << indexOfFirstTwoSiteOp << " "
                        << indexOfSecondTwoSiteOp << std::endl;
            else
                infoout << "No" << std::endl;
            infoout << "Observables threshold: " << observableThreshold
                    << std::endl;
        };
    }
    else
    {
        std::cerr << "Couldn't open output files." << std::endl;
        exit(EXIT_FAILURE);
    };
    infoout.close();
    stepData data; // this struct will contain most of the important parameters
    data.compBlock = data.beforeCompBlock = NULL;
    for(int trial = 1; trial <= nTrials; trial++)
    {
        clock_t startTrial = clock();
        std::cout << "Trial " << trial << ":" << std::endl;
        std::ofstream fileout("Output/Trial_" + std::to_string(trial));
        fileout << "Trial " << trial << ":\n" << std::endl;
        
        // read in parameters that vary over trials:
        int lSys;                           // system length
        filein >> lSys;
        std::vector<double> couplingConstants(nCouplingConstants);
        for(int i = 0; i < nCouplingConstants; i++)
            filein >> couplingConstants[i];
        double k;              // wave number of rotation of interstitial spins
        int rangeOfObservables, // number of sites at which to measure observables
            nSweeps;                        // number of sweeps to be performed
        filein >> k >> rangeOfObservables >> data.mMax >> nSweeps;
        if(rangeOfObservables == -1)
            rangeOfObservables = lSys;
        std::vector<double> groundStateErrorTolerances(nSweeps + 1);
        for(int sweep = 0; sweep <= nSweeps; sweep++)
            filein >> groundStateErrorTolerances[sweep];
        
        fileout << "System length: " << lSys << "\nCoupling constants:";
        for(double couplingConstant : couplingConstants)
            fileout << " " << couplingConstant;
        fileout << "\nWave number: " << k << "\nMaximum bond dimension: "
                << data.mMax << "\nNumber of sweeps: " << nSweeps
                << "\nLanczos tolerances:";
        for(double groundStateErrorTolerance : groundStateErrorTolerances)
            fileout << " " << groundStateErrorTolerance;
        fileout << std::endl << std::endl;
        data.ham.setParams(couplingConstants, lSys, k);
        int skips = 0,
            runningKeptStates = d * d;
        for(; runningKeptStates <= data.mMax; skips++)
            runningKeptStates *= d;  // find how many edge sites can be skipped
        bool oddSize = lSys % 2;
        int lSFinal,                        // final length of the system block
            lEFinal;                   // final length of the environment block
        if(oddSize)
        {
            lSFinal = (lSys - 1)/2;
            lEFinal = (lSys - 3)/2;
        }
        else
            lSFinal = lEFinal = lSys / 2 - 1;
        bool completeED = false;
        if(skips + 1 >= lSFinal)
        {
            if(skips + 1 == lSFinal && runningKeptStates == data.mMax * d)
            {
                std::cout << "Note: the maximum bond dimension is large enough "
                          << "to perform exact diagonalization." << std::endl;
                completeED = true;
            }
            else
            {
                std::cout << "Error: the maximum bond dimension is larger than "
                          << "required for exact diagonalization." << std::endl;
                continue;
            };
        };
        std::vector<TheBlock> westBlocks(lSys - 2 - skips),
                              eastBlocks(lSys - 2 - skips);
             // initialize system - the last block is only used for odd-size ED
        TheBlock* eastBlocksStart = eastBlocks.data();
        westBlocks.front() = TheBlock(data.ham, true);
        eastBlocks.front() = TheBlock(data.ham, false);
                                         // initialize the edge one-site blocks
        std::cout << "Performing iDMRG..." << std::endl;
            // note: this iDMRG code assumes parity symmetry of the Hamiltonian
        data.sweepingEast = true;
        data.exactDiag = true;
        data.compBlock = eastBlocksStart;
        data.infiniteStage = true;
        data.lancTolerance = groundStateErrorTolerances.front()
                             * groundStateErrorTolerances.front() / 2;
        rmMatrixX_t psiGround;                    // seed for Lanczos algorithm
        double cumulativeTruncationError = 0.;
        for(int site = 0; site < skips; site++, data.compBlock++) // initial ED
            westBlocks[site + 1]
                = westBlocks[site].nextBlock(data, psiGround,
                                             cumulativeTruncationError,
                                             eastBlocksStart + site + 1);
        data.exactDiag = completeED;
        for(int site = skips, end = lEFinal - 1; site < end; site++,
                                                             data.compBlock++)
                                                                       // iDMRG
        {
            psiGround = randomSeed(westBlocks[site], eastBlocks[site]);
            westBlocks[site + 1]
                = westBlocks[site].nextBlock(data, psiGround,
                                             cumulativeTruncationError,
                                             eastBlocksStart + site + 1);
        };
        if(oddSize)
        {
            data.compBlock = eastBlocksStart + (lSFinal - 2);
            psiGround = randomSeed(westBlocks[lSFinal - 2],
                                   eastBlocks[lSFinal - 2]);
            westBlocks[lSFinal - 1]
                = westBlocks[lSFinal - 2].nextBlock(data, psiGround,
                                                    cumulativeTruncationError);
        };
        std::cout << "iDMRG average truncation error: "
                  << cumulativeTruncationError / (lSys - 2)
                  << std::endl;       // handles both even and odd system sizes
        if(completeED || nSweeps == 0)
            psiGround = randomSeed(westBlocks[lSFinal - 1],
                                   eastBlocks[lEFinal - 1]);
        else
        {
            std::cout << "Performing fDMRG..." << std::endl;
            data.infiniteStage = false;
            int endSweep = lSys - 4 - skips;              // last site of sweep
            psiGround = randomSeed(westBlocks[lSFinal - 1],
                                   eastBlocks[lEFinal - 1]);
            for(int sweep = 1; sweep <= nSweeps; sweep++)
                                                    // perform the fDMRG sweeps
            {
                data.compBlock = eastBlocksStart + (lEFinal - 1);
                data.lancTolerance = groundStateErrorTolerances[sweep]
                                     * groundStateErrorTolerances[sweep] / 2;
                data.beforeCompBlock = data.compBlock - 1;
                cumulativeTruncationError = 0.;
                for(int site = lSFinal - 1; site < endSweep;
                    site++, data.compBlock--, data.beforeCompBlock--)
                    westBlocks[site + 1]
                        = westBlocks[site].nextBlock(data, psiGround,
                                                     cumulativeTruncationError);
                reflectPredictedPsi(psiGround, westBlocks[endSweep],
                                    eastBlocks[skips]);
                               // reflect the system to reverse sweep direction
                data.sweepingEast = false;
                data.compBlock = &westBlocks[endSweep];
                data.beforeCompBlock = data.compBlock - 1;
                for(int site = skips; site < endSweep;
                    site++, data.compBlock--, data.beforeCompBlock--)
                    eastBlocks[site + 1]
                        = eastBlocks[site].nextBlock(data, psiGround,
                                                     cumulativeTruncationError);
                reflectPredictedPsi(psiGround, eastBlocks[endSweep],
                                    westBlocks[skips]);
                data.sweepingEast = true;
                data.compBlock = eastBlocksStart + endSweep;
                data.beforeCompBlock = data.compBlock - 1;
                for(int site = skips, end = lSFinal - 1; site < end;
                    site++, data.compBlock--, data.beforeCompBlock--)
                    westBlocks[site + 1]
                        = westBlocks[site].nextBlock(data, psiGround,
                                                     cumulativeTruncationError);
                std::cout << "Sweep " << sweep
                          << " complete. Average truncation error: "
                          << cumulativeTruncationError / (4 * lSys - 8)
                          << std::endl;
            };
        };
        data.compBlock = eastBlocksStart + (lEFinal - 1);
        data.infiniteStage = false;
        FinalSuperblock hSuperFinal
            = westBlocks[lSFinal - 1].createHSuperFinal(data, psiGround, skips);
                                               // calculate ground-state energy
        fileout << "Ground state energy density = "
                << hSuperFinal.gsEnergy / lSys << std::endl << std::endl;
        if(calcObservables)
        {
            std::cout << "Calculating observables..." << std::endl;
            VectorXd oneSiteVals;
            MatrixXd twoSiteVals;
            if(calcOneSiteExpValues)   // calculate one-site expectation values
                oneSiteVals = oneSiteExpValues(oneSiteOp, rangeOfObservables,
                                               lSys, hSuperFinal, westBlocks,
                                               eastBlocks, fileout);
            if(calcTwoSiteExpValues)   // calculate two-site expectation values
                twoSiteVals = twoSiteExpValues(firstTwoSiteOp, secondTwoSiteOp,
                                               rangeOfObservables, lSys,
                                               hSuperFinal, westBlocks,
                                               eastBlocks, fileout);
            if(calcOneSiteExpValues && calcTwoSiteExpValues)
            {
                MatrixXd connectedCorrFunc
                    = twoSiteVals - oneSiteVals * oneSiteVals.transpose();
                for(int i = 0, end = rangeOfObservables * rangeOfObservables;
                    i < end; i++)
                    if(std::abs(connectedCorrFunc(i)) < observableThreshold)
                        connectedCorrFunc(i) = 0.;
                fileout << "Connected correlation function:\n"
                        << connectedCorrFunc << std::endl << std::endl;
            };
        };
        std::cout << std::endl;
        clock_t stopTrial = clock();
        fileout << "Elapsed time: "
                << float(stopTrial - startTrial)/CLOCKS_PER_SEC << " s"
                << std::endl;
        fileout.close();
    };
    filein.close();
    
    clock_t stop = clock();
    std::cout << "Done. Elapsed time: " << float(stop - start)/CLOCKS_PER_SEC
              << " s" << std::endl;

    return 0;
}
