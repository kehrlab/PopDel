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
#include "popdel_view/popdel_view.h"
#include "popdel_view/popdel_view_parameter_parsing.h"
#include "popdel_call/translocation_popdel_call.h"
#include "popdel_call/detect_junction_popdel_call.h"

using namespace seqan;

typedef Iterator<String<String<Call> >, Rooted>::Type TBufferIt;               // Iterator over the buffer
typedef Iterator<String<bool> >::Type TWWCIt;

// ==========================================================================
// Function processSegment()
// ==========================================================================
// Controlls the functions for the window-wise genotyping of one chromosome.
void processSegment(ChromosomeProfile & chromosomeProfile,
                    TranslocationProfile & translocProfile,
                    String<Call> & deletionCalls,
                    Tuple<String<JunctionCall>, 4> & junctionCalls,
                    Windows & windows,
                    VcfOutputBundle & vcfOutput,
                    PopDelCallParameters & params,
                    const TRGs & rgs)
{
    // Check if start and endProfiles are in sync and sync them otherwise.
    // If we reached the end of the startProfiles during the last call of processSegment, the profiles are
    // lacking behing by one window. This fixes it
    if (chromosomeProfile.profilesAtEnd)
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));

    chromosomeProfile.initializeActiveReads();
    chromosomeProfile.initializeActiveLoad(rgs);
    // translocProfile.setWindow(chromosomeProfile.chrom, chromosomeProfile.currentPos);
    while (!chromosomeProfile.profilesAtEnd)
    {
        genotype_deletion_window(deletionCalls, chromosomeProfile, rgs, params);
        if (!params.delOnly)
        {
            bool junctionCallMade = detect_junction_window(junctionCalls[0],        // rfCalls
                                                        junctionCalls[1],        // ffCalls
                                                        junctionCalls[2],        // rrCalls
                                                        junctionCalls[3],        // translocCalls
                                                        chromosomeProfile,
                                                        translocProfile,
                                                        rgs,
                                                        params);
            assignCoverageChanges(junctionCalls[0],
                                junctionCalls[1],
                                junctionCalls[2],
                                chromosomeProfile,
                                rgs,
                                junctionCallMade);
            windows.currentWindow.i2 += params.windowShift;
        }
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));
        // translocProfile.setWindow(chromosomeProfile.chrom, chromosomeProfile.currentPos);
    }
    // if (!empty(junctionCalls))
    // {
    //     printJunctionCalls(junctionCalls);
    // }
    bool jncFound = false;
    if (!params.delOnly)
        jncFound = mergeJunctionWindows(junctionCalls, params.meanStddev);

    bool delFound = mergeDeletionWindows(deletionCalls, params.meanStddev, params.minRelWinCover, params.outputFailed);
    // if (jncFound)
    // {
    //     printJunctionCalls(junctionCalls);
    // }
    if (delFound || jncFound)
    {
        vcfOutput.appendContigName(windows.currentWindow.i1);
        writeRegenotypedCalls(vcfOutput, deletionCalls, junctionCalls, params);
        vcfOutput.flush();
    }
}
struct WindowWise
{};
// Overload for window-wise outpu.
void processSegment(ChromosomeProfile & chromosomeProfile,
                    TranslocationProfile & translocProfile,
                    String<Call>  & calls,
                    Tuple<String<JunctionCall>, 4> & junctionCalls,
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
    bool junctionCallMade = false;

    while (!chromosomeProfile.profilesAtEnd)                                        // Now move over the chromosome window-by-window.
    {
        genotype_deletion_window(calls, chromosomeProfile, rgs, params);
        if  (!params.delOnly)
        {
            junctionCallMade |= detect_junction_window(junctionCalls[0],
                                                       junctionCalls[1],
                                                       junctionCalls[2],
                                                       junctionCalls[3],
                                                       chromosomeProfile,
                                                       translocProfile,
                                                       rgs,
                                                       params);
        }
        windows.currentWindow.i2 += params.windowShift;
        chromosomeProfile.profilesAtEnd = !(chromosomeProfile.nextWindow(params.windowShift));
        // translocProfile.setWindow(chromosomeProfile.chrom, chromosomeProfile.currentPos);
    }
    if (!empty(calls))
         vcfOutput.appendContigName(windows.currentWindow.i1);

    writeRegenotypedCalls(vcfOutput, calls, junctionCalls, params);
    vcfOutput.flush();

    // if (junctionCallMade)    // For testing purposes
    //     printJunctionCalls(junctionCalls);
}
// ==========================================================================
// Function markRoiAsCompleteForFile()
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

    String<CharString> cNames;
    getContigNames(cNames, bwi.infile);
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
               length(cNames),
               params.mergeRG);

    // Open the output file and a compression stream.
    std::ofstream out(toCString(params.outfile), std::ios::out | std::ios::binary);
    if (!out.good())
        SEQAN_THROW(FileOpenError(toCString(params.outfile)));

    // Initialize file offsets (by chromosome) for profile index.
    unsigned indexSize = 0;
    String<int32_t> cLengths;
    getContigLengths(cLengths, bwi.infile);


    bwi.blacklist.loadFromBED(params.translocationBlacklist, cNames, getMaxHistReadLen(params.histograms));

    String<String<uint64_t> > indexFields;              // Prepare index structure
    resize(indexFields, length(cLengths), Exact());
    for (unsigned i = 0; i < length(cLengths); ++i)
    {
        resize(indexFields[i], cLengths[i]/params.indexRegionSize + 1, 0);
        indexSize += length(indexFields[i]);
    }
    String<uint64_t> translocationIndexFields;         // Prepare translocation index structure
    resize(translocationIndexFields, length(cLengths), 0, Exact());

    // Write uncompressed header to the output file.
    writeProfileHeader(out,
                       params.sampleName,
                       params.indexRegionSize,
                       indexSize,
                       length(translocationIndexFields),
                       params.readGroups,
                       params.histograms,
                       cNames,
                       cLengths,
                       params.mergeRG);
    bool bwiNext = true;
    bwi.qualReq.sameChrom = false;  // After creating the header we also allow translocated reads.
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
            for (unsigned rg = 0; rg < length((*bwi).records); ++rg)
                std::sort(begin((*bwi).records[rg]), end((*bwi).records[rg]));

            // Write function subtracts mean insert size from all records in window.
            writeWindow(zipper, *bwi);
            bwiNext = goNext(bwi);
        }
        while (bwiNext && regionChrom == (*bwi).chrom && regionPos == (*bwi).beginPos / (__int64)params.indexRegionSize);
        zipper.zflush();
    }
    writeTranslocationBlock(out, translocationIndexFields, bwi.translocationBuffers);
    uint64_t maxOffset = out.tellp();
    // Write the profile index into the header.
    writeIndexIntoHeader(out, indexFields, maxOffset);
    writeTranslocationIndexIntoHeader(out,
                                      translocationIndexFields,
                                      maxOffset);
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
    String<Call> bufferedDeletions;        //create and resize output buffer.
    reserve(bufferedDeletions, params.windowBuffer / 10);
    Tuple<String<JunctionCall>, 4> bufferedJunctions;
    for (unsigned i = 0; i < 4; ++i)
        reserve(bufferedJunctions[i], params.windowBuffer / 10);


    // Objects for keeping track of position on genome.
    Windows windows;
    windows.nextWindow = windows.currentWindow = getFirstWindowCoordinate(bufferedWindows, finishedROIs, params);
    String<Pair<CharString, uint32_t> > nextCandidateWindows;    // Holds the last read positions of each sample.
    resize(nextCandidateWindows, params.sampleNum, windows.currentWindow, Exact());
    unsigned nextReadPos = maxValue<unsigned>();

    // Create vcf output object.
    VcfOutputBundle vcfOutput(params.outfile, params.sampleNames, params.contigNames[0], params.contigLengths[0]);
    vcfOutput.appendContigName(windows.currentWindow.i1);

    // Create the object for managing all the read pairs from all samples and set it to the first position.
    ChromosomeProfile chromosomeProfile(length(params.readGroups), params.maxLoad, params.windowBuffer);
    TranslocationProfile translocProfile(length(params.readGroups));
    chromosomeProfile.currentPos = windows.currentWindow.i2;
    chromosomeProfile.resetTo(chromosomeProfile.currentPos);

    while (params.fileCount != 0)           // filecount is reduced by one for every finished file.
    {
        if (!goNextRegion(windows, params))
            return 0;

        for (unsigned s = 0; s < params.sampleNum; ++s)             // TODO: Move this loop into seperate function
        {
            if (finishedFiles[s] || finishedROIs[s])               // Don't try to read finished files.
                continue;

            std::ifstream inStream(toCString(params.inputFiles[s]), std::ios::in | std::ios::binary);
            if (!inStream.good())
                SEQAN_THROW(FileOpenError(toCString(params.inputFiles[s])));

            if (goNextRoi(inStream, s, params.nextRoi, finishedROIs, params))      // Try to go to the next Roi.
            {
                if (checkAndSwitch(chromosomeProfile, params.rgs[s], nextCandidateWindows[s].i2))
                {
                    unsigned contigIdx = params.representativeContigs ? 0 : s;
                    unsigned segmentCode = readSegment(chromosomeProfile,        // 0 - Segment finished
                                                       inStream,                 // 1 - ROI or Chrom not yet reached
                                                       params.inputFiles[s],     // 2 - End of ROI or chromosome
                                                       params.rgs[s],            // 3 - End of file
                                                       params.histograms,
                                                       params.contigNames[contigIdx],
                                                       nextCandidateWindows[s],
                                                       *params.nextRoi,
                                                       bufferedWindows[s]);
                    // if (inStream.eof())
                    //         inStream.clear();

                    // uint64_t translocIndexBeginPos = translocationIndexBeginPos(params.numRegions[s]);
                    // jumpToTranslocRegion(inStream,
                    //                      params.contigNames[contigIdx],
                    //                      params.nextRoi->seqName,
                    //                      translocIndexBeginPos);
                    // readTranslocationSegment(translocProfile,
                    //                          inStream,
                    //                          params.inputFiles[s],
                    //                          params.rgs[s],
                    //                          params.contigNames[contigIdx],
                    //                          *params.nextRoi);

                    if (!processSegmentCode(nextReadPos,
                                            finishedROIs,
                                            finishedFiles,
                                            params,
                                            windows,
                                            s,
                                            nextCandidateWindows,
                                            segmentCode))
                    {
                        continue;
                    }
                    inStream.close();
                }
            }
            else    // ROI is done for this sample.
            {
                markRoiAsCompleteForFile(finishedROIs, finishedFiles, params.fileCount, s);
            }
        }
        if (!params.windowWiseOutput)
            processSegment(chromosomeProfile,
                           translocProfile,
                           bufferedDeletions,
                           bufferedJunctions,
                           windows,
                           vcfOutput,
                           params,
                           params.rgs);
        else
            processSegment(chromosomeProfile,
                           translocProfile,
                           bufferedDeletions,
                           bufferedJunctions,
                           windows,
                           vcfOutput,
                           params,
                           params.rgs,
                           WindowWise());

        clear(bufferedDeletions);
        for (unsigned i = 0; i < 4; ++i)
            clear(bufferedJunctions[i]);

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
    CharString sampleName;
    uint16_t profileVersion;
    String<CharString> readGroups;
    String<Histogram> histograms;
    String<CharString> contigNames;
    String<int32_t> contigLengths;
    unsigned numRegions;
    unsigned indexRegionSize;
    readProfileHeader(in,
                      profileVersion,
                      params.infile,
                      sampleName,
                      readGroups,
                      histograms,
                      contigNames,
                      contigLengths,
                      numRegions,
                      indexRegionSize);
    const unsigned numReadGroups = length(readGroups);
    // Write the header to std::cout.
    if (params.writeHeader || params.writeOnlyHeader)
        writeProfileHeader(std::cout, profileVersion, sampleName, readGroups, contigNames, contigLengths);

    if (params.writeHistograms)
    {
        printHistograms(std::cout, readGroups, histograms);
        return 0;
    }

    if (params.writeOnlyHeader)
        return 0;

    if (params.region.seqName == "")
    {
        if (params.translocationsOnly)
        {
            uint64_t translocIndexBeginPos = translocationIndexBeginPos(numRegions);
            jumpToTranslocBlock(in, translocIndexBeginPos);
            printTranslocationWindows(in, params, contigNames, numReadGroups);
        }
        else
        {
            printWindows(in, params, contigNames, numReadGroups);
            if (params.viewTranslocations)
            {
                if (in.eof())
                    in.clear();

                std::cout << "#Read pairs mapping to different contigs:" << std::endl;
                uint64_t translocIndexBeginPos = translocationIndexBeginPos(numRegions);
                jumpToTranslocBlock(in, translocIndexBeginPos);
                printTranslocationWindows(in, params, contigNames, numReadGroups);
            }
        }
    }
    else
    {
        fillInvalidPositions(params.region, contigNames, contigLengths);
        if (params.translocationsOnly)
        {
            uint64_t translocIndexBeginPos = translocationIndexBeginPos(numRegions);
            jumpToTranslocRegion(in, contigNames, params.region.seqName, translocIndexBeginPos);
            printTranslocationWindowsInRegion(in, params, contigNames, numReadGroups);
        }
        else
        {
            // Read file offset from index and move the stream there.
            jumpToRegion(in, contigNames, contigLengths, indexRegionSize, params.region);
            printWindowsInRegion(in, params, contigNames, numReadGroups);
            if (params.viewTranslocations)
            {
                if (in.eof())
                    in.clear();

                std::cout << "#Read pairs mapping to different contigs:" << std::endl;
                uint64_t translocIndexBeginPos = translocationIndexBeginPos(numRegions);
                jumpToTranslocRegion(in, contigNames, params.region.seqName, translocIndexBeginPos);
                printTranslocationWindowsInRegion(in, params, contigNames, numReadGroups);
            }
        }
    }
    return 0;
}

#endif /*WORKFLOW_POPDEL_H_*/
