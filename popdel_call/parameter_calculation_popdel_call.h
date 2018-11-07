#ifndef PARAMETER_CALCULATION_POPDEL_CALL_H_
#define PARAMETER_CALCULATION_POPDEL_CALL_H_

#include "utils_popdel.h"
#include "parameter_parsing_popdel_call.h"
#include "../insert_histogram_popdel.h"
#include "load_profile_popdel_call.h"       //for initializeRois

using namespace seqan;

// -----------------------------------------------------------------------------
// Function calculateMinimalLikelihoodRatio()
// -----------------------------------------------------------------------------
// Calculates the default threshold for a significant call.
// Uses the 95% quantile of the chi squared distribution with 1 degree of freedom.
inline void calculateMinimalLikelihoodRatio(double & minLikelihoodRatio, double prior) //TODO: Check DF
{
    SEQAN_ASSERT_GT(prior, 0);
    SEQAN_ASSERT_LT(prior, 1);
    double const c = 6.6349;                        //99% quantile of chi squared distribution, df=1. Right tail only.
    minLikelihoodRatio = (c / 2.0) - log(prior / (1 - prior));
    SEQAN_ASSERT_GT(minLikelihoodRatio, 0);
}
// -----------------------------------------------------------------------------
// Function setMinimumDeletionLengths()
// -----------------------------------------------------------------------------
// Set the minimum length of deletions to the value of minDelLen[0].
// If not specified, set it to the standard deviation of the insert sizes for each read group.
inline void setMinInitDelLengths(String<unsigned> & minDelLen,
                                 const String<Histogram> & histograms,
                                 double factor = 4) // one-sided 99.9% quantile
{
    SEQAN_ASSERT_GT(factor, 0);
    std::ostringstream msg;
    if (length(minDelLen) == 1)
    {
        resize(minDelLen, length(histograms), minDelLen[0]);
        msg << "The minimum initial deletion length has been set to " << minDelLen[0] << " for all read groups.";
    }
    else
    {
        resize(minDelLen, length(histograms));
        unsigned minLen = maxValue<unsigned>();
        unsigned maxLen = 0;
        for (unsigned i = 0; i < length(histograms); ++i)
        {
            minDelLen[i] = _round(factor * histograms[i].stddev);
            if (minLen > minDelLen[i])
                minLen = minDelLen[i];
            if (maxLen < minDelLen[i])
                maxLen = minDelLen[i];
        }
        msg << "Minimum initial deletion lengths have been set to " << factor << " * standard deviations of insert size histograms [" << minLen 
        << ".." << maxLen << "].";
    }
    printStatus(msg);
}
// =======================================================================================
// Function setMinDelLengths()
// =======================================================================================
void setMinDelLength(PopDelCallParameters & params, const double & factor = 1.0)
{
    SEQAN_ASSERT_GT(factor, 0);
    if (params.minLen == maxValue<unsigned>())  // Only do calculation when minLen not set by user.
        params.minLen = _round(factor * calculatePercentile(params.minInitDelLengths, 0.95));
}
// =======================================================================================
// Function setAvgStddev()
// =======================================================================================
// Calculate the average of the histograms' standard deviations and set the variable in params.
void setMeanStddev(PopDelCallParameters & params)
{
    SEQAN_ASSERT_GEQ(length(params.histograms), 1u);
    double sum = 0;
    for(Iterator<String<Histogram> >::Type it = begin(params.histograms); it != end(params.histograms); ++it)
        sum += it->stddev;
    params.meanStddev = sum / length(params.histograms);
}
// =======================================================================================
// Function loadAndCalculateParameters()
// =======================================================================================
// Load all insert size histograms and perform the parameter estimation (or set them as specified).
inline void loadAndCalculateParameters(PopDelCallParameters & params)
{
    // Load the insert size histograms and profile headers.
    loadInsertSizeHistograms(params.histograms,
                             params.readGroups,
                             params.rgs,
                             params.inputFiles,
                             params.contigNames,
                             params.contigLengths,
                             params.indexRegionSizes,
                             params.smoothing);
    std::ostringstream msg;
    msg << "Loaded insert size histograms for " << length(params.histograms) << " read groups.";
    printStatus(msg);
    if (params.windowShift == 0)
        params.windowShift = params.windowSize;
    setMeanStddev(params);
    setMinInitDelLengths(params.minInitDelLengths, params.histograms);
    setMinDelLength(params);
    calculateMinimalLikelihoodRatio(params.minimumLikelihoodRatio, params.prior);
    if (!empty(params.roiFile))
    {
        initializeRois(params.allRois, params.nextRoi, params.roiList, params.contigNames[0], params.roiFile);
        checkRois(params.allRois, params.contigNames[0]);
    }
    else if (!empty(params.roiList))
    {
        initializeRois(params.allRois, params.nextRoi, params.roiList, params.contigNames[0]);
        checkRois(params.allRois, params.contigNames[0]);
    }
    else
    {   // If there are no defined ROIs: Take the contigs of the first sample.
        std::vector<std::string> roiList;
        append(roiList, params.contigNames[0]);
        initializeRois(params.allRois, params.nextRoi, roiList, params.contigNames[0]);
    }
    msg.str("");
    msg << "Calculated minimum likelihood ratio as " << params.minimumLikelihoodRatio << " from the prior probability " 
        << params.prior << " using the 95%-quantile of a chi-squared distribution with df=1.";
    printStatus(msg);
}

#endif /* PARAMETER_CALCULATION_POPDEL_CALL_H_ */
