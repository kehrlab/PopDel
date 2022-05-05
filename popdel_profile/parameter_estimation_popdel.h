#ifndef PARAMETER_ESTIMATION_POPDEL_H_
#define PARAMETER_ESTIMATION_POPDEL_H_

#include <fstream>
#include <seqan/bam_io.h>

#include "../utils_popdel.h"
#include "../histogram_popdel.h"
#include "profile_parameter_parsing_popdel.h"

using namespace seqan;


void printParams(const CharString & rg, const unsigned & n, const Histogram & hist)
{
    std::cout << "All read groups merged into \'" << rg << "\':";
    std::cout << "  SampledReadPairs=" << n;
    std::cout << "  Median=" << hist.median;
    std::cout << "  Mean=" << hist.mean;
    std::cout << "  StdDev=" << hist.stddev;
    std::cout << "  AvgCoverage=" << hist.coverage;
    std::cout << "  ReadLength=" << hist.readLength;
    std::cout << std::endl;
}
// ---------------------------------------------------------------------------------------
// Function printReadGroupParams()
// ---------------------------------------------------------------------------------------
// Print information about each read group to standard output.
void printReadGroupParams(const PopDelProfileParameters & params)
{
    typedef PopDelProfileParameters::TReadGroups::const_iterator TReadGroupIter;
    TReadGroupIter it = params.readGroups.begin();
    TReadGroupIter itEnd = params.readGroups.end();
    if (params.mergeRG && params.readGroups.size() > 1)
        printParams(it->first, params.sampleSize[it->second], params.histograms[it->second]);
    else
        for (; it != itEnd; ++it)
            printParams(it->first, params.sampleSize[it->second], params.histograms[it->second]);
}
// ---------------------------------------------------------------------------------------
// Function printWindowParams()
// ---------------------------------------------------------------------------------------
// Print window-size and shift per iteration to standard output.
void printWindowParameters(const PopDelProfileParameters & params)
{
    std::ostringstream msg;
    msg << "Window parameters:";
    msg << "  WindowSize=" << params.windowSize;
    msg << "  WindowShift=" << params.windowShift;
    printStatus(msg);
}
// ---------------------------------------------------------------------------------------
// Function getReadLength()
// ---------------------------------------------------------------------------------------
// Return the read length infered from a BamAlignmentRecord.
// Corrects for hard-clipped bases by adding them to the length of the sequence
inline unsigned getReadLength(const BamAlignmentRecord & record)
{
    SEQAN_ASSERT_GT(record._n_cigar, 0u);
    unsigned clip = 0;
    if (record.cigar[0].operation == 'H')
        clip += record.cigar[0].count;

    if (record.cigar[record._n_cigar - 1].operation == 'H')
        clip += record.cigar[record._n_cigar - 1].count;

    return record._l_qseq + clip;
}
// Read all records in the specified region of infile, check and add them to the set of insert size histograms.
inline void addToHistograms(String<Histogram> & histograms,
                            const GenomicRegion & interval,
                            HtsFile & infile,
                            String<unsigned> & sampleSize,
                            BamQualReq & qualReq,
                            const std::map<CharString, unsigned> & readGroups)
{
    bool mergeRG = (length(histograms) == 1u && readGroups.size() > 1);
    if (!setRegion(infile, toCString(interval.seqName), interval.beginPos, interval.endPos))
        return;

    BamAlignmentRecord record;
    while (seqFreeReadRegion(record, infile))
    {
        if (record.rID == interval.rID && record.beginPos < interval.beginPos)      //Iterate till specified region.
            continue;
        if (record.rID != interval.rID || record.beginPos > interval.endPos)        //Break if region has been passed.
            break;
        if (meetsRequirements(record, qualReq) && getPairOrientation(record.flag) == Orientation::FR)
        {
            unsigned rg = getReadGroup(record.tags, readGroups, mergeRG);
            if (histograms[rg].readLength == 0)
                histograms[rg].readLength = getReadLength(record);

            unsigned insertSize = std::min(std::abs(record.tLen),
                                           static_cast<__int32>(length(histograms[rg].values) - 1));
            histograms[rg].values[insertSize] += 1.0;
            ++sampleSize[rg];
        }
    }
}
// ---------------------------------------------------------------------------------------
// Function getHistograms()
// ---------------------------------------------------------------------------------------
typedef Iterator<const String<GenomicRegion> >::Type TItvIter;
// Call addToHistograms for each specified interval.
// Return the iterator to the first unused interval (or to the end of the container).
inline TItvIter getHistograms (String<Histogram> & histograms,
                               HtsFile & infile,
                               const String<GenomicRegion> & intervals,
                               PopDelProfileParameters & params)
{

    Histogram hist;
    resize(hist.values, 20000, 0.0);
    if (params.mergeRG)
    {
        resize(histograms, 1, hist, Exact());
        resize(params.sampleSize, 1, 0, Exact());
    }
    else
    {
        resize(histograms, params.readGroups.size(), hist, Exact());
        resize(params.sampleSize, length(params.readGroups), 0, Exact());
    }
    TItvIter itv = begin(intervals);
    TItvIter itvEnd = end(intervals);

    bool enough = false;
    while (itv != itvEnd && !enough)                               //Call _addToHistograms for each interval.
    {
        addToHistograms(histograms,
                        *itv,
                        infile,
                        params.sampleSize,
                        params.qualReq,
                        params.readGroups);
        enough = true;
        for (Iterator<const String<unsigned> >::Type it = begin(params.sampleSize); it != end(params.sampleSize); ++it)
        {
            if (*it < params.minSampling)
            {
                enough = false;
                break;
            }
        }
        ++itv;
    }
    if (!enough)
    {
        std::ostringstream msg;
        msg << "[ERROR] [PopDel] Not enough reads in sampling regions for read group(s) ";
        bool p = false;
        for (PopDelProfileParameters::TReadGroups::const_iterator it = params.readGroups.begin();
             it != params.readGroups.end(); ++it)
        {
            if (params.sampleSize[it->second] < params.minSampling)
            {
                if (p)
                    msg << ", ";
                else
                    p = true;

                msg << "\'";
                msg <<  it->first;
                msg << "\'";
            }
        }
        msg << ". Please use different sampling intervals (option '-i') or use the option '-n'"
            << " to reduce the number of required read pairs for each histogram of each read group.";
        SEQAN_THROW(IOError(toCString(msg.str())));
    }
    else
    {
        std::ostringstream msg;
        if (params.mergeRG && params.readGroups.size() > 1)
            msg << "Determined inner distance histogram of merged read groups using the interval(s): ";
        else
            msg << "Determined inner distance histograms for all read groups using the interval(s): ";
        bool p = false;
        for (TItvIter it = begin(intervals); it != itv; ++it)
        {
            if (p)
                msg << ", ";
            else
                p = true;

            msg << it->seqName << ":" << it->beginPos + 1 << "-" << it->endPos;
        }
        printStatus(msg);
    }
    return (itv);
}
// ---------------------------------------------------------------------------------------
// Function totalIntervalLength()
// ---------------------------------------------------------------------------------------
//Calculate the total length of all intervals by adding their individual lengths.
inline __uint32 totalIntervalLength(const String<GenomicRegion> & intervals, const TItvIter & itvEnd)
{
    typedef Iterator<const String<GenomicRegion> >::Type TItvIter;
    __uint32 totalItvLen = 0;
    TItvIter itv = begin(intervals);
    while (itv != itvEnd)
    {
        totalItvLen += ((*itv).endPos - (*itv).beginPos) - 1;
        ++itv;
    }
    return totalItvLen;
}
// ---------------------------------------------------------------------------------------
// Function calculateCoverage()
// ---------------------------------------------------------------------------------------
// Calculate the (average) coverage for a single read group
inline void calculateCoverage(Histogram & histogram, __uint32 totalIntervalLength)
{
    typedef Iterator<const String<double>, Rooted>::Type THistIter;
    double counts = 0;
    THistIter it = begin(histogram.values, Rooted());
    THistIter itEnd = end(histogram.values);
    while (it != itEnd)
    {
        counts += *it;
        ++it;
    }
//     std::cout << "sampledReads: " << counts << " totalItvLength: "
//               << totalIntervalLength << " readLength: " << histogram.readLength << std::endl;
    histogram.coverage = (counts / totalIntervalLength) * histogram.readLength * 2;
}
// =======================================================================================
// Function checkCoverage()
// =======================================================================================
// Gives a warning of the coverage of a read group is below minCov
inline void checkCoverage(const PopDelProfileParameters & params,
                          const unsigned minCov = 5)
{
    typedef PopDelProfileParameters::TReadGroups::const_iterator TReadGroupIter;
    TReadGroupIter itEnd = params.readGroups.end();
    bool p = false;
    for (TReadGroupIter it = params.readGroups.begin(); it != itEnd; ++it)
    {
        if (params.histograms[it->second].coverage < minCov)
        {
            p = true;
            std::ostringstream msg;
            if (params.mergeRG  && params.readGroups.size() > 1)
            {
                msg << "[WARNING] Coverage of merged read groups is below " << minCov
                    << " across the sampled intervals (see above).";
                printStatus(msg);
                break;
            }
            else
            {
                msg << "[WARNING] Coverage of read group\'" <<  it->first << "\' is below " << minCov
                    << " across the sampled intervals (see above).";
                     printStatus(msg);
            }
        }
    }
    if (p)
        std::cout << std::endl;
}
// =======================================================================================
// Function calculateParamters()
// =======================================================================================
// Sets the parameters for each read group.
void calculateParameters(HtsFile & infile, PopDelProfileParameters & params)
{
    if (params.outfile == "*BAM/CRAM-FILE*.profile")
    {
        params.outfile = params.bamfile;
        append(params.outfile, ".profile");
    }
    getReadGroups(params.readGroups, infile, params.mergeRG); // Get the read groups and chromosome names/lengths from the file header.
    if (params.sampleName == "")
        if (!getSampleName(params.sampleName, infile))
            params.sampleName = getSampleNameFromPath(params.bamfile);
    std::ostringstream msg;
    msg << "Found " << length(params.readGroups) << " read groups in input alignment file \'" << params.bamfile << "\'.";
    if (params.mergeRG && params.readGroups.size() > 1)
        msg << " All read groups will be merged into read group \'" << params.readGroups.begin()->first << "\'.";

    printStatus(msg);
    if (length(params.rois) == 0)                                      // Estimate the genomic regions.
    {
        unsigned totalLength = getWholeGenomeIntervals(params.rois, infile);
        std::ostringstream msg;
        msg << "Reference consists of " << length(params.rois)
            << " sequences with a total length of " << totalLength << " base pairs.";
        printStatus(msg);
    }
    String<GenomicRegion> intervals;
    readIntervals(intervals, params.intervalFile, infile, params.rois, params.referenceVersion);
    mergeOverlappingIntervals(intervals);                                // Handle the regular intervals for sampling
    setRIDs(intervals, infile);
    std::stable_sort(begin(intervals), end(intervals), &lowerRIDGenomicRegion);

    TItvIter itvEnd = getHistograms (params.histograms, infile, intervals, params);// Calculate hist for each RG...
                                                                        //...and get the read length at the same time.
    unsigned totalItvLen = totalIntervalLength(intervals, itvEnd);
    for (unsigned i = 0; i < length(params.histograms); ++i)            // For all histograms...
    {                                                   // Calculate mean, stddev, and coverage if not user-specified.
        calculateCoverage(params.histograms[i], totalItvLen);
        calculateHistMetrics(params.histograms[i]);
    }
    printStatus("Histogram statistics (median, mean, stddev) calculated.");
    //writeHistograms(params.outfile, params.readGroups, params.histograms);
    if (params.windowShift == 0)             // Set the window shift to window size if it was not specified by the user.
        params.windowShift = params.windowSize;
}
#endif /* PARAMETER_ESTIMATION_POPDEL_H_ */
