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
        resize(minDelLen, length(histograms), minDelLen[0], Exact());
        msg << "The minimum initial deletion length has been set to " << minDelLen[0] << " for all samples and readgroups.";
    }
    else
    {
        resize(minDelLen, length(histograms), Exact());
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
    std::ostringstream msg;
    if (params.minLen == maxValue<unsigned>())  // Only do calculation when minLen not set by user.
        params.minLen = _round(factor * calculatePercentile(params.minInitDelLengths, 0.95));
    msg << "Minimum final deletion length has been set to " << params.minLen << ".";
    printStatus(msg);
}
// =======================================================================================
// Function setAvgStddev()
// =======================================================================================
// Calculate the average of the histograms' standard deviations and set the variable in params.
void setMeanStddev(PopDelCallParameters & params)
{
    SEQAN_ASSERT_GEQ(length(params.histograms), 1u);
    double sum = 0;
    for(Iterator<String<const Histogram> >::Type it = begin(params.histograms); it != end(params.histograms); ++it)
        sum += it->stddev;
    params.meanStddev = sum / length(params.histograms);
}
// =======================================================================================
// Function loadMaxLoad()
// =======================================================================================
inline bool loadMaxLoad(PopDelCallParameters & params)
{
    if (params.defaultMaxLoad == 0u)
        params.defaultMaxLoad = maxValue<unsigned>();
    resize(params.maxLoad, params.readGroups.size(), params.defaultMaxLoad);
    if (empty(params.maxLoadFile))
        return true;
    std::ifstream file(toCString(params.maxLoadFile));
    std::ostringstream msg;
    if (!file.is_open())
    {
        msg << "[PopDel] Could not open coverage file \'" << params.maxLoadFile << "\' for reading.";
        SEQAN_THROW(IOError(toCString(msg.str())));
        return false;
    }
    std::string readGroup;
    std::string load;
    unsigned len = 0;
    std::map<CharString, unsigned> loadMap;
    while (file >> readGroup)
    {
        if (file >> load)
        {
            loadMap[readGroup] = std::stoul(load);
            ++len;
        }
        else
        {
            msg << "[PopDel] Could not read maximum coverage for read group \'" << readGroup << "\'.";
            SEQAN_THROW(IOError(toCString(msg.str())));
            return false;
        }
    }
    // For every RG, look if there is a specified maxCov in the loadMap. Keep the default value for the rest.
    unsigned countAssigned = 0;
    String<CharString> unassignedRG;
    for (std::map<CharString, unsigned>::const_iterator it = params.readGroups.begin();
         it != params.readGroups.end();
         ++it)
    {
        std::map<CharString, unsigned>::const_iterator lIt = loadMap.find(it->first);
        if (lIt != loadMap.end())
        {
            if (lIt->second == 0u)
                params.maxLoad[it->second] = maxValue<unsigned>();
            else
                params.maxLoad[it->second] = lIt->second;
            ++countAssigned;
        }
        else
        {
            appendValue(unassignedRG, it->first);
        }
    }
    msg << "Loaded maximum active coverage for " << len << " read groups. (" << countAssigned << " assigned)";
    printStatus(msg);
    if (!empty(unassignedRG))
    {
        msg.str("");
        msg << "WARNING: No maximum active coverage value specified for the following read groups: ";
        for (Iterator<const String<CharString> >::Type it = begin(unassignedRG); it != end(unassignedRG); ++it)
        {
            if (it != begin(unassignedRG))
                msg << ", ";
            msg << *it;
        }
        msg << "\nAssuming the default value of " << params.defaultMaxLoad << " for the listed read group(s).";
        printStatus(msg);
    }
    return true;
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
                             params.smoothing,
                             params.modRgByFileName,
                             params.representativeContigs,
                             params.pseudoCountFraction);
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
    loadMaxLoad(params);
    msg.str("");
    msg << "Calculated minimum log-likelihood ratio as " << params.minimumLikelihoodRatio << " from the prior probability "
        << params.prior << " using the 99%-quantile of a chi-squared distribution with df=1.";
    printStatus(msg);
}

#endif /* PARAMETER_CALCULATION_POPDEL_CALL_H_ */
