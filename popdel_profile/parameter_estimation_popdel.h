#ifndef PARAMETER_ESTIMATION_POPDEL_H_
#define PARAMETER_ESTIMATION_POPDEL_H_

#include <fstream>
#include <seqan/bam_io.h>

#include "../utils_popdel.h"
#include "../insert_histogram_popdel.h"
#include "profile_parameter_parsing_popdel.h"

using namespace seqan;

// ---------------------------------------------------------------------------------------
// Function printReadGroupParams()
// ---------------------------------------------------------------------------------------
// Print information about each read group to standard output.
void printReadGroupParams(const PopDelProfileParameters & params)
{
    typedef PopDelProfileParameters::TReadGroups::const_iterator TReadGroupIter;
    TReadGroupIter itEnd = params.readGroups.end();
    if (params.mergeRG && params.readGroups.size() > 1)
    {
        TReadGroupIter it = params.readGroups.begin();
        std::cout << "All read groups merged into \'" << it->first << "\':";
        std::cout << "  SampledReadPairs=" << params.sampleSize[it->second];
        std::cout << "  Median=" << params.histograms[it->second].median;
        std::cout << "  Mean=" << params.histograms[it->second].mean;
        std::cout << "  StdDev=" << params.histograms[it->second].stddev;
        std::cout << "  AvgCoverage=" << params.histograms[it->second].coverage;
        std::cout << "  ReadLength=" << params.histograms[it->second].readLength;
        std::cout << std::endl;
    }
    else
    {
        for (TReadGroupIter it = params.readGroups.begin(); it != itEnd; ++it)
        {
            std::cout << "Read group \'" << it->first << "\':";
            std::cout << "  SampledReadPairs=" << params.sampleSize[it->second];
            std::cout << "  Median=" << params.histograms[it->second].median;
            std::cout << "  Mean=" << params.histograms[it->second].mean;
            std::cout << "  StdDev=" << params.histograms[it->second].stddev;
            std::cout << "  AvgCoverage=" << params.histograms[it->second].coverage;
            std::cout << "  ReadLength=" << params.histograms[it->second].readLength;
            std::cout << std::endl;
        }
    }
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
// Function addToInsertSizeHistograms()
// ---------------------------------------------------------------------------------------
// Read all records in the specified region of infile, check and add them to the set of insert size histograms.
inline void addToInsertSizeHistograms(String<Histogram> & histograms,
                                      const GenomicRegion & interval,
                                      BamFileIn & infile,
                                      const BamIndex<Bai> & bai,
                                      String<unsigned> & sampleSize,
                                      BamQualReq & qualReq,
                                      const std::map<CharString, unsigned> & readGroups)
{
    bool mergeRG = (length(histograms) == 1u && readGroups.size() > 1);
    bool alis;
    jumpToRegion(infile, alis, interval.rID, interval.beginPos, interval.endPos, bai);
    if (!alis)
        return;
    BamAlignmentRecord record;
    while (!atEnd(infile))
    {
        try
        {
            readRecord(record, infile);
        }
        catch (...)
        {
            SEQAN_THROW(IOError("Could no read record in BAM-File."));
        }
        if (record.rID == interval.rID && record.beginPos < interval.beginPos)      //Iterate till specified region.
            continue;
        if (record.rID != interval.rID || record.beginPos > interval.endPos)        //Break if region has been passed.
            break;
        if (!meetsRequirements(record, qualReq) || record.beginPos > record.pNext)  //Check requirements of record.
            continue;
        unsigned rg = getReadGroup(record.tags, readGroups, mergeRG);
        unsigned insertSize = std::min(std::abs(record.tLen), static_cast<__int32>(length(histograms[rg].values) - 1));
        histograms[rg].values[insertSize] += 1.0;
        ++sampleSize[rg];
        if (histograms[rg].readLength == 0)
        {
            if (record.cigar[0].operation == 'H')
            {
                unsigned clip = record.cigar[0].count;
                if (record.cigar[length(record.cigar) - 1].operation != 'H')
                {
                    clip += record.cigar[length(record.cigar) - 1].count;
                }
                histograms[rg].readLength = length(record.seq) + clip;
            }
            else
            {
                histograms[rg].readLength = length(record.seq);
            }
        }
    }
}
// ---------------------------------------------------------------------------------------
// Function getInsertSizeHistograms()
// ---------------------------------------------------------------------------------------
typedef Iterator<const String<GenomicRegion> >::Type TItvIter;
// Call addToInsertSizeHistograms for each specified interval.
// Return the iterator to the first unused interval (or to the end of the container).
inline TItvIter getInsertSizeHistograms(String<Histogram> & histograms,
                                        BamFileIn & infile,
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
    BamIndex<Bai> bai;
    loadBai(bai, params.bamfile);
    TItvIter itv = begin(intervals);
    TItvIter itvEnd = end(intervals);

    bool enough = false;
    while (itv != itvEnd && !enough)                               //Call _addToInsertSizeHistograms for each interval.
    {
        addToInsertSizeHistograms(histograms,
                                  *itv, infile,
                                  bai,
                                  params.sampleSize,
                                  params.qualReq,
                                  params.readGroups);
        enough = true;
        for (Iterator<String<unsigned> >::Type it = begin(params.sampleSize); it != end(params.sampleSize); ++it)
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
            msg << "Determined insert size histogram of merged read groups using the interval(s): ";
        else
            msg << "Determined insert size histograms for all read groups using the interval(s): ";
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
// ---------------------------------------------------------------------------------------
// Function writeHistograms()
// ---------------------------------------------------------------------------------------
// Write all histograms to output file, one line per read goup.
// Contents per line: RG min max median stddev readLength counts1 counts2 ... countsN
void writeHistograms(const CharString & outfile,
                     const std::map<CharString, unsigned> & readGroups,
                     const String<Histogram> & histograms)
{
    typedef std::map<CharString, unsigned>::const_iterator TReadGroupIter;
    typedef Iterator<const String<double> >::Type THistIter;
    CharString histogramFile = outfile;
    histogramFile += ".hist";
    // Open output file.
    std::ofstream histfile(toCString(histogramFile));
    if (!histfile.is_open())
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not open file \'" << histogramFile << "\' for writing.";
        SEQAN_THROW(IOError(toCString(msg.str())));
    }
    // Write histogram on one line per read group.
    for(TReadGroupIter rgIt = readGroups.begin(); rgIt != readGroups.end(); ++rgIt)
    {
        const Histogram & hist = histograms[rgIt->second];
        unsigned minValue = getHistLeftBorder(hist);
        unsigned maxValue = getHistRightBorder(hist);
        histfile << rgIt->first << "\t" << minValue << "\t" << maxValue;
        histfile << "\t" << hist.median << "\t" << hist.stddev << "\t" << hist.readLength;
        THistIter it = begin(hist.values);
        THistIter itEnd = begin(hist.values);
        it += minValue;
        itEnd += maxValue;
        while (it != itEnd)
        {
            histfile << "\t" << *it;
            ++it;
        }
        histfile << std::endl;
    }
    std::ostringstream msg;
    msg << "Histograms written to \'" << histogramFile << "\'.";
    printStatus(msg);
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
void calculateParameters(PopDelProfileParameters & params)
{
    if (params.outfile == "*BAM-FILE*.profile")
    {
        params.outfile = params.bamfile;
        append(params.outfile, ".profile");
    }
    BamFileIn infile;
    open(infile, toCString(params.bamfile));
    BamHeader header;
    readHeader(header, infile);
    getReadGroups(params.readGroups, header, params.mergeRG); // Get the read groups and chromosome names/lengths from bam file header.
    std::ostringstream msg;
    msg << "Found " << length(params.readGroups) << " read groups in input bam file \'" << params.bamfile << "\'.";
    if (params.mergeRG && params.readGroups.size() > 1)
        msg << " All read groups will be merged into read group \'" << params.readGroups.begin()->first << "\'.";
    printStatus(msg);
    if (length(params.rois) == 0)                                      // Estimate the genomic regions.
    {
        unsigned totalLength = getWholeGenomeIntervals(params.rois, header);
        std::ostringstream msg;
        msg << "Reference consists of " << length(params.rois)
            << " sequences with a total length of " << totalLength << " base pairs.";
        printStatus(msg);
    }
      // Read list of intervals to use for parameter calculation.
    String<GenomicRegion> intervals;
    readIntervals(intervals, params.intervalFile, header, params.rois);
    mergeOverlappingIntervals(intervals);                                // Handle the regular intervals for sampling
    setRIDs(intervals, infile);
    TItvIter itvEnd = getInsertSizeHistograms(params.histograms, infile, intervals, params);// Calculate insert size hist for each RG...
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
