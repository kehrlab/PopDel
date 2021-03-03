#ifndef WORKFLOW_POPDEL_H_
#define WORKFLOW_POPDEL_H_

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

#include "parse_popdel.h"
#include "utils_popdel.h"
#include "popdel_profile/profile_parameter_parsing_popdel.h"
#include "popdel_profile/parameter_estimation_popdel.h"
#include "popdel_profile/bam_window_iterator_popdel.h"
#include "popdel_call/parameter_calculation_popdel_call.h"
#include "popdel_call/genotype_deletion_popdel_call.h"
#include "popdel_call/load_profile_popdel_call.h"
#include "popdel_call/vcfout_popdel_call.h"
#include "popdel_view_parameter_parsing.h"

using namespace seqan;

typedef Iterator<String<String<Call> >, Rooted>::Type TBufferIt;               // Iterator over the buffer
typedef Iterator<String<bool> >::Type TWWCIt;

// ==========================================================================
// Function processSegment()
// ==========================================================================
// Controlls the functions for the window-wise genotyping of one chromosome.
void processSegment(ChromosomeProfile & chromosomeProfile,
                    String<Call> & calls,
                    Windows & windows,
                    VcfOutputBundle & vcfOutput,
                    PopDelCallParameters & params,
                    const TRGs & rgs)
{
    // Check if start and endProfiles are in sync and sync them otherwise.
    // If we reached the end of the startProfiles during the last call of processSegment, the profiles are
    // lacking behing by one window. This fixes it.
    if (chromosomeProfile.profilesAtEnd)
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));

    chromosomeProfile.initializeActiveReads();
    while (!chromosomeProfile.profilesAtEnd)
    {
        genotype_deletion_window(calls, chromosomeProfile, rgs, params);
        windows.currentWindow.i2 += params.windowShift;
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));
    }
    if(unifyCalls(calls, params.meanStddev, params.minRelWinCover, params.outputFailed))
    {
        vcfOutput.appendContigName(windows.currentWindow.i1);
        writeRegenotypedCalls(vcfOutput, calls, params);
        vcfOutput.flush();
    }
}
struct WindowWise
{};
// Overload for window-wise outpu.
void processSegment(ChromosomeProfile & chromosomeProfile,
                    String<Call>  & calls,
                    Windows & windows,
                    VcfOutputBundle & vcfOutput,
                    PopDelCallParameters & params,
                    const TRGs & rgs,
                    WindowWise ww)
{
    (void) ww;
    // Check if start and endProfiles are in sync and sync them otherwise.
    // If we reached the end of the startProfiles during the last call of processChromosome, the profiles are
    // lacking behing by one window. This fixes it.
    if (chromosomeProfile.profilesAtEnd)
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));

    chromosomeProfile.initializeActiveReads();

    while (!chromosomeProfile.profilesAtEnd)                                        // Now move over the chromosome window-by-window.
    {
        genotype_deletion_window(calls, chromosomeProfile, rgs, params);
        windows.currentWindow.i2 += params.windowShift;
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));
    }
    if (!empty(calls))
         vcfOutput.appendContigName(windows.currentWindow.i1);

    writeRegenotypedCalls(vcfOutput, calls, params);
    vcfOutput.flush();
}
// ==========================================================================
// Function processSegment()
// ==========================================================================
// Mark the ROI as completely processed for the given file.
// If the file has been completely processed mark it and decrease the counter for unfinished files.
inline void markRoiAsCompleteForFile(String<bool> & finishedROIs,
                                     String<bool> & finishedFiles,
                                     unsigned & fileCount,
                                     const unsigned & sampleNum)
{
    finishedROIs[sampleNum] = true;
    if (!finishedFiles[sampleNum])
    {
        finishedFiles[sampleNum] = true;
        SEQAN_ASSERT_GT(fileCount, 0u);
        --fileCount;
    }
}
// =======================================================================================
// Function finalizeRoi()
// =======================================================================================
// Perfom all neccessary steps for finishing the processing of the current ROI and advancing to the next one.
// Return false if the last ROI has been processed, true otherwise.
inline bool finalizeRoi(ChromosomeProfile & chromosomeProfile,
                        PopDelCallParameters & params,
                        unsigned & nextReadPos,
                        Windows & windows,
                        VcfOutputBundle & vcfOutput,
                        String<bool> & finishedROIs)
{
    if (!advanceRoiIterator(params.nextRoi, params.allRois))
    {   //Last ROI has been processed.
        std::ostringstream msg;
        msg << "Output written to \'" << params.outfile << "\'.";
        printStatus(msg);
        return false;
    }
    else
    {
        chromosomeProfile.fullReset();
        nextReadPos = maxValue<unsigned>();
        windows.currentWindow.i1 = params.nextRoi->seqName;
        windows.currentWindow.i2 = 0;
        windows.nextChromWindow.i1 = "";
        windows.nextChromWindow.i2 = maxValue<unsigned>();
        vcfOutput.unlockContigName();
        for (Iterator<String<bool> >::Type it = begin(finishedROIs); it != end(finishedROIs); ++it)
            *it = false;

        return true;
    }
}
// =======================================================================================
// Function popdel_profile()
// =======================================================================================
// Wrapper for calling all neccessary functions for performing the profile creation.
int popdel_profile(int argc, char const ** argv)
{
    // Parse the command line and set all user-specified parameters.
    //   Print a header and the command to std::out.
    PopDelProfileParameters params;
    int res = checkParser(parseCommandLine(params, argc, argv));
    if (res >= 0)
        return res;

    // Calculate parameters for the input data:
    //   Get read groups from bam file, genome intervals, calculate insert size
    BamWindowIterator bwi(params.bamfile, params.referenceFile);
    //   distributions, estimate window size, etc.
    calculateParameters(bwi.infile, params);
    printStatus("Finished parameter estimation.\n");
    printReadGroupParams(params);
    std::cout << std::endl;
    checkCoverage(params);

    // Define an iterator over the whole bam file or the specified intervals and iterate window by window:
    //   Calculate standardized read pair distance histograms for each genomic window
    //   and write the corresponding insert size profile to the output file.

    initialize(bwi,
               params.readGroups,
               params.rois,
               params.qualReq,
               params.maxDeletionSize,
               params.windowSize,
               params.windowShift,
               params.histograms,
               params.mergeRG);

    // Open the output file and a compression stream.
    std::ofstream out(toCString(params.outfile), std::ios::out | std::ios::binary);
    if (!out.good())
        SEQAN_THROW(FileOpenError(toCString(params.outfile)));

    // Initialize file offsets (by chromosome) for profile index.
    unsigned indexSize = 0;
    String<int32_t> cLengths;
    getContigLengths(cLengths, bwi.infile);

    String<CharString> cNames;
    getContigNames(cNames, bwi.infile);

    String<String<uint64_t> > indexFields;
    resize(indexFields, length(cLengths));
    for (unsigned i = 0; i < length(cLengths); ++i)
    {
        resize(indexFields[i], cLengths[i]/params.indexRegionSize + 1, 0);
        indexSize += length(indexFields[i]);
    }

    // Write uncompressed header to the output file.
    writeProfileHeader(out,
                       params.indexRegionSize,
                       indexSize,
                       params.readGroups,
                       params.histograms,
                       cNames,
                       cLengths,
                       params.mergeRG);

    bool bwiNext = true;
    while (bwiNext)
    {
        // Register file offset of region for index.
        __int32 regionChrom = (*bwi).chrom;
        __int32 regionPos = (*bwi).beginPos / params.indexRegionSize;
        indexFields[regionChrom][regionPos] = out.tellp();

        // Write windows of the next region as a compressed block.
        zlib_stream::zip_ostream zipper(out);
        do
        {
            // Sort the insert sizes for each read group before writing them to output.
            for (unsigned rg = 0; rg < length((*bwi).insertSizes); ++rg)
                std::sort(begin((*bwi).insertSizes[rg]), end((*bwi).insertSizes[rg]));

            // Write function subtracts mean insert size from all records in window.
            writeWindow(zipper, *bwi, params.histograms);
            bwiNext = goNext(bwi);
        }
        while (bwiNext && regionChrom == (*bwi).chrom && regionPos == (*bwi).beginPos / (__int64)params.indexRegionSize);
        zipper.zflush();
    }

    // Write the profile index into the header.
    writeIndexIntoHeader(out, indexFields, out.tellp());

    std::ostringstream msg;
    msg << "Output written to \'" << params.outfile << "\'.";
    printStatus(msg);
    return 0;
}
// =======================================================================================
// Function popdel_call()
// =======================================================================================
// Wrapper for calling all neccessary functions for performing the deletion calling.
int popdel_call(int argc, char const ** argv)
{
    // Parse the command line.
    PopDelCallParameters params;
    int res = checkParser(parseCommandLine(params, argc, argv));
    if (res >= 0)
        return res;

    //Calculate all parameters
    loadAndCalculateParameters(params);
    printStatus("Finished parameter calculation.\n");
    std::ostringstream msg;
    msg << "Initialized insert size profiles from " << params.sampleNum << " input files.";
    printStatus(msg);

    // Create objects for tracking the progress and buffering.
    String<bool> finishedFiles;             // True for every sample completely done.
    resize(finishedFiles, params.sampleNum, false, Exact());
    String<bool> finishedROIs;        // True for every sample that finished the current ROI.
    resize(finishedROIs, params.sampleNum, false, Exact());
    String<String<Window> > bufferedWindows;
    resize(bufferedWindows, params.sampleNum, Exact());
    String<Call> bufferedCalls;        //create and resize output buffer.
    reserve(bufferedCalls, params.windowBuffer / 10);

    // Objects for keeping track of position on genome.
    Windows windows;
    windows.nextWindow = windows.currentWindow = getFirstWindowCoordinate(bufferedWindows, finishedROIs, params);
    String<Pair<CharString, uint32_t> > nextCandidateWindows;    // Holds the last read positions of each sample.
    resize(nextCandidateWindows, params.sampleNum, windows.currentWindow, Exact());
    unsigned nextReadPos = maxValue<unsigned>();

    // Create vcf output object.
    VcfOutputBundle vcfOutput(params.outfile, params.inputFiles, params.contigNames[0], params.contigLengths[0]);
    vcfOutput.appendContigName(windows.currentWindow.i1);

    // Create the object for managing all the read pairs from all samples and set it to the first position.
    ChromosomeProfile chromosomeProfile(length(params.readGroups), params.maxLoad, params.windowBuffer);
    chromosomeProfile.currentPos = windows.currentWindow.i2;
    chromosomeProfile.resetTo(chromosomeProfile.currentPos);

    while (params.fileCount != 0)           // filecount is reduced by one for every finished file.
    {
        if (!goNextRegion(windows, params))
            return 0;

        for (unsigned i = 0; i < params.sampleNum; ++i)
        {
            if (finishedFiles[i] || finishedROIs[i])               // Don't try to read finished files.
                continue;

            std::ifstream inStream(toCString(params.inputFiles[i]), std::ios::in | std::ios::binary);
            if (!inStream.good())
                SEQAN_THROW(FileOpenError(toCString(params.inputFiles[i])));

            if (goNextRoi(inStream, i, params.nextRoi, finishedROIs, params))      // Try to go to the next Roi.
            {
                if (checkAndSwitch(chromosomeProfile, params.rgs[i], nextCandidateWindows[i].i2))
                {
                    unsigned contigIdx = params.representativeContigs?0:i;
                    unsigned segmentCode = readSegment(chromosomeProfile,        // 0 - Segment finished
                                                       inStream,                 // 1 - ROI or Chrom not yet reached
                                                       params.inputFiles[i],     // 2 - End of ROI or chromosome
                                                       params.rgs[i],            // 3 - End of file
                                                       params.histograms,
                                                       params.contigNames[contigIdx],
                                                       nextCandidateWindows[i],
                                                       *params.nextRoi,
                                                       bufferedWindows[i]);
                    inStream.close();
                    if (!processSegmentCode(nextReadPos,
                                            finishedROIs,
                                            finishedFiles,
                                            params,
                                            windows,
                                            i,
                                            nextCandidateWindows,
                                            segmentCode))
                        continue;
                }
            }
            else    // ROI is done for this sample.
            {
                markRoiAsCompleteForFile(finishedROIs, finishedFiles, params.fileCount, i);
            }
        }
        if (!params.windowWiseOutput)
            processSegment(chromosomeProfile, bufferedCalls, windows, vcfOutput, params, params.rgs);
        else
            processSegment(chromosomeProfile, bufferedCalls, windows, vcfOutput, params, params.rgs, WindowWise());

        clear(bufferedCalls);
        if (allTrue(finishedROIs))  // ROI done for all samples.
        {
            if(!finalizeRoi(chromosomeProfile, params, nextReadPos, windows, vcfOutput, finishedROIs))
                return 0;           // Last ROI has been processed.

            getFirstWindowOnNextROI(windows.currentWindow, bufferedWindows, finishedROIs, params);
            chromosomeProfile.currentPos = windows.currentWindow.i2;
            chromosomeProfile.resetTo(chromosomeProfile.currentPos);
            nextReadPos = maxValue<unsigned>();
            for (auto it = begin(nextCandidateWindows); it != end(nextCandidateWindows); ++it)
                *it = windows.currentWindow;
        }
        else
        {
            windows.currentWindow.i2 = nextReadPos;
            nextReadPos = maxValue<unsigned>();
        }
    }
    msg.str("");
    msg << "Output written to \'" << params.outfile << "\'.";
    printStatus(msg);
    return 0;
}

// =======================================================================================
// Function popdel_view()
// =======================================================================================
int popdel_view(int argc, char const ** argv)
{
    PopDelViewParameters params;
    int res = checkParser(parseCommandLine(params, argc, argv));
    if (res >= 0)
        return res;

    // Open input file and decompression stream.
    std::ifstream in(toCString(params.infile), std::ios::in | std::ios::binary);
    if (!in.good())
        SEQAN_THROW(FileOpenError(toCString(params.infile)));

    // Read the header.
    String<CharString> readGroups;
    String<Histogram> histograms;
    String<CharString> contigNames;
    String<int32_t> contigLengths;
    unsigned indexRegionSize;
    readProfileHeader(in, params.infile, readGroups, histograms, contigNames, contigLengths, indexRegionSize);
    // Write the header to std::cout.
    if (params.writeHeader || params.writeOnlyHeader)
        writeProfileHeader(std::cout, readGroups, contigNames, contigLengths);

    if (params.writeHistograms)
    {
        printHistograms(std::cout, readGroups, histograms);
        return 0;
    }

    if (params.writeOnlyHeader)
        return 0;

    Window window;
    if (params.region.seqName == "")
    {
        // Read all windows and write them to std::cout.
        zlib_stream::zip_istream unzipper(in);
        while(readWindow(unzipper, window, length(readGroups)))
        {
            writeWindow(std::cout, window, contigNames);
        }
    }
    else
    {
        // Read file offset from index and move the stream there.
        fillInvalidPositions(params.region, contigNames, contigLengths);
        jumpToRegion(in, contigNames, contigLengths, indexRegionSize, params.region);

        // Write all windows in the region to std::cout.
        zlib_stream::zip_istream unzipper(in);
        while(readWindow(unzipper, window, length(readGroups)))
        {
            if (contigNames[window.chrom] != params.region.seqName ||
               (contigNames[window.chrom] == params.region.seqName &&
                          window.beginPos >= params.region.endPos))
                break;
            if (params.region.beginPos > window.beginPos)
                continue;
            writeWindow(std::cout, window, contigNames);
        }
    }

    return 0;
}

#endif /*WORKFLOW_POPDEL_H_*/
