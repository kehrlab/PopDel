#ifndef LOAD_PROFILE_POPDEL_CALL_H_
#define LOAD_PROFILE_POPDEL_CALL_H_

#include <seqan/sequence.h>

#include "../utils_popdel.h"
#include "profile_structure_popdel_call.h"
#include "parameter_parsing_popdel_call.h"
#include "../popdel_profile/window_podel.h"

using namespace seqan;

typedef Iterator<const String<GenomicRegion>, Rooted>::Type TRegionIter;
// ======================================================================================
// Class Windows
// ======================================================================================
// Class for holding/managing all information on the windows.
struct Windows
{
    Pair<CharString, __uint32> currentWindow;           // The current Window
    Pair<CharString, __uint32> nextWindow;              // The potential next window
    Pair<CharString, __uint32> nextChromWindow;         // The potential next window on the next chromosome.

    Windows() :
    currentWindow("", 0),
    nextWindow("", 0),
    nextChromWindow("", maxValue<__uint32>())
    {}

    inline bool nextChromosome()
    {
        if (currentWindow.i1 != nextChromWindow.i1)
        {
            currentWindow = nextChromWindow;
            nextWindow = nextChromWindow;
            nextChromWindow.i1 = "";
            return true;
        }
        else
        {
            nextWindow = currentWindow;
            return false;
        }
    }
};
// =======================================================================================
// Function checkRois()
// =======================================================================================
// Check if all ROIs are valid, i.e. they appear in the list of contig names.
// Throw an error if one of the ROIs is not contained in contigNames.
// Otherwise return true;
bool checkRois(const String<GenomicRegion> & allRois, const String<CharString> & contigNames)
{
    Iterator<const String<GenomicRegion>, Standard>::Type roiIt = begin(allRois);
    Iterator<const String<GenomicRegion>, Standard>::Type roiEnd = end(allRois);
    Iterator<const String<CharString>,  Standard>::Type namesEnd =end(contigNames);
    while (roiIt != roiEnd)
    {
        Iterator<const String<CharString>,  Standard>::Type namesIt = begin(contigNames);
        while (namesIt != namesEnd)
        {
            if (*namesIt == roiIt->seqName)
                break;
            else
                ++namesIt;
        }
        if (namesIt == namesEnd)
        {
            std::ostringstream msg;
            msg << "[PopDel] Invalid chromosome name \'" << roiIt->seqName <<
                   "\' not part of reference sequence names in BAM header. Terminating.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
            return false;
        }
        ++roiIt;
    }
    return true;
}
// Reorders the ROIS in allRois s.t. the order of the sequence names matches the one in desiredOrder.
inline void recreateOrder(String<GenomicRegion> & allRois, const String<String<char> > & desiredOrder)
{
    String<Pair<unsigned>> order;
    resize(order, length(allRois), Exact());
    for (unsigned i = 0; i < length(allRois); ++i)
    {
        for (unsigned j = 0; j < length(desiredOrder); ++j)
        {
            if (allRois[i].seqName == desiredOrder[j])
            {
                order[i] = Pair<unsigned>(j, i);
                break;
            }
        }
    }
    std::sort(begin(order), end(order));
    String<GenomicRegion> reOrdered;
    reserve(reOrdered, length(allRois), Exact());
    for (Iterator<String<Pair<unsigned> > >::Type it = begin(order); it != end(order); ++it)
        append(reOrdered, allRois[it->i2]);

    move(allRois, reOrdered);
}
// =======================================================================================
// Function initializeRois()
// =======================================================================================
// Prepare the list of regions of interest per chromosome.
// Return true on success, false otherwise or if the set is empty.
inline bool initializeRois(String<GenomicRegion> & allRois,
                           Iterator<String<GenomicRegion> >::Type & nextRoi,
                           const std::vector<std::string> & intervalStrings,
                           const String<String<char> >  & contigNames) //Needed for correct ordering.
{
    typedef std::vector<std::string>::const_iterator TIter;
    typedef std::map<CharString, String<GenomicRegion> >::iterator TRoiIter;
    std::map<CharString, String<GenomicRegion> > roiMap;

    // Sort the ROIs by chromsome
    for (TIter it = intervalStrings.begin(); it != intervalStrings.end(); ++it)
    {
        GenomicRegion itv;
        if(!parseGenomicRegion(itv, *it))
        {
            std::ostringstream msg;
            msg << "[PopDel] Error while parsing genomic region \'" << *it << "\'. Terminating.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
            return false;
        }

        if (itv.beginPos < 0)
            itv.beginPos = 0;
        if (itv.endPos < 0)
            itv.endPos = maxValue<__int32>();

        appendValue(roiMap[itv.seqName], itv);
    }
    unsigned c = 0; // Get total number of ROIs.
    for (TRoiIter rIt = roiMap.begin(); rIt != roiMap.end(); ++rIt)
        c += length(rIt->second);

    if (c == 0)
        return false;

    resize(allRois, c);
    unsigned i = 0;
    for (TRoiIter rIt = roiMap.begin(); rIt != roiMap.end(); ++rIt)
    {
        mergeOverlappingIntervals(rIt->second);
        for (unsigned j = 0; j < length(rIt->second); ++j)
        {
            allRois[i] = (rIt->second[j]);
            ++i;
        }
    }
    recreateOrder(allRois, contigNames);
    nextRoi = begin(allRois);
    return true;
}
//Overload for loading ROIs only from file or adding them to those loaded from command line.
inline bool initializeRois(String<GenomicRegion> & allRois,
                           Iterator<String<GenomicRegion> >::Type & nextRoi,
                           std::vector<std::string> & intervalStrings,
                           const String<String<char> >  & contigNames,
                           const CharString & filename)
{
    std::ifstream file(toCString(filename));
    if (!file.is_open())
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not open interval file \'" << filename << "\' for reading.";
        SEQAN_THROW(IOError(toCString(msg.str())));
        return false;
    }
    std::string word;
    unsigned len = 0;
    while (file >> word)
    {
        intervalStrings.push_back(word);
        ++len;
    }
    std::ostringstream msg;
    msg << "Finished reading " << len << " intervals from file \'" << filename << "\'.";
    printStatus(msg);
    return initializeRois(allRois,nextRoi, intervalStrings, contigNames);
}
// =======================================================================================
// Function peekNextRoiIterator()
// =======================================================================================
// Return false if the last ROI has been processed, true otherwise.
inline bool peekNextRoiIterator(const Iterator<String<GenomicRegion> >::Type & nextRoi,
                                const String<GenomicRegion> & allRois)
{
    SEQAN_ASSERT_NEQ(nextRoi, end(allRois));
    if (nextRoi + 1 == end(allRois))
        return false;
    else
        return true;
}
// =======================================================================================
// Function advanceRoiIterator()
// =======================================================================================
// Advance nextRoi.
// Return false if the last ROI has been processed, true otherwise.
inline bool advanceRoiIterator(Iterator<String<GenomicRegion> >::Type & nextRoi,
                               const String<GenomicRegion> & allRois)
{
    SEQAN_ASSERT_NEQ(nextRoi, end(allRois));
    ++nextRoi;
    if (nextRoi == end(allRois))
        return false;
    else
        return true;
}
// ======================================================================================
// Function adaptRegions()
// ======================================================================================
// Adapt the next ROI s.t. it matches the leftmost next window candidate or adapt the next candiate window to match
// the ROI depending on what comes first.
// Return false if the next ROI lies on another chromosome or if the next candidate window comes after the ROI.
// Return true otherwise.
inline bool adaptRegions(Pair<CharString, unsigned> & nextWindowCandidate, GenomicRegion & nextRoi, unsigned w)
{
    if (nextWindowCandidate.i1 != nextRoi.seqName)
        return false;
    else if (static_cast<int>(nextWindowCandidate.i2) >= nextRoi.endPos)
        return false;
    else if (static_cast<int>(nextWindowCandidate.i2 + w - 1)== nextRoi.beginPos) // Nothing to do.
        return true;
    else if (static_cast<int>(nextWindowCandidate.i2 + w - 1) > nextRoi.beginPos)
        nextRoi.beginPos = nextWindowCandidate.i2;              // Adjust the left end of the ROI.
    else
        nextWindowCandidate.i2 += ((nextRoi.beginPos - nextWindowCandidate.i2) / w) * w;  // Move the candidate to the beginning of the ROI.

    return true;
}
// =======================================================================================
// Function goNextRoi()
// =======================================================================================
// Advance nextRoi if all samples are done with the current ROI and jump to its region.
// Return false if the last ROI has been processed, true otherwise.
inline bool goNextRoi(std::ifstream & file,
                      const unsigned & fileNum,
                      Iterator<String<GenomicRegion> >::Type & nextRoi,
                      String<bool> & finishedROIs,
                      const PopDelCallParameters & params)
{
    if (allTrue(finishedROIs))
    {
        if (!advanceRoiIterator(nextRoi, params.allRois))
            return false;
        else
            for (Iterator<String<bool> >::Type it = begin(finishedROIs); it != end(finishedROIs); ++it)
                *it = false;
    }
    unsigned idx = params.representativeContigs?0:fileNum;
    jumpToRegion(file,
                 params.contigNames[idx],
                 params.contigLengths[idx],
                 params.indexRegionSizes[idx],
                 *nextRoi);
    return true;
}
// ======================================================================================
// Function goNextRegion()
// ======================================================================================
// Advance to the next region of interest or set the current ROI to continue where we left of.
// Return false if no more ROIs are left, true otherwise.
inline bool goNextRegion(Windows & windows, PopDelCallParameters & params)
{
    while (!adaptRegions(windows.currentWindow, *params.nextRoi, params.windowSize))
    {   // No sample contains the contig of the next ROI. Skip ROIs until we find the next valid contig.
        while (advanceRoiIterator(params.nextRoi, params.allRois))
        {
            if (params.nextRoi->seqName  == windows.currentWindow.i1 &&
                params.nextRoi->beginPos >= static_cast<int>(windows.currentWindow.i2))
                break;
        }
        if (params.nextRoi == end(params.allRois))
        {   // No matching ROIs left. Terminate
            std::ostringstream msg;
            msg << "Output written to \'" << params.outfile << "\'.";
            printStatus(msg);
            return false;
        }
    }
    return true;
}
// ======================================================================================
// Function getFirstWindowCoordinate()
// ======================================================================================
// Return the first (=smallest) window position from all input streams.
Pair<CharString, __uint32> getFirstWindowCoordinate (String<String<Window> >& convertedWindows,
                                                     String<bool> & finishedROIs,
                                                     PopDelCallParameters & params)
{
    Pair<CharString, __uint32> minCoord("", maxValue<__uint32>());
    Pair<CharString, __uint32> currentCoord;
    std::string str;
    Window window;
    while (true)
    {
        for (unsigned i = 0; i < params.fileCount; ++i)
        {
            unsigned contigIdx = params.representativeContigs?0:i;
            String<Window> & currentSampleConvertedWindows = convertedWindows[i];
            std::ifstream in(toCString(params.inputFiles[i]), std::ios::in | std::ios::binary);
            if (!in.good())
                SEQAN_THROW(FileOpenError(toCString(params.inputFiles[i])));

            //jump to the first region of interest.
            goNextRoi(in, i, params.nextRoi, finishedROIs, params);

            readWindow(in, window, length(params.rgs[i]), params.uncompressedIn);

            if (isEmpty(window))
                continue;

            //Convert the 256bp window to a string of 30bp windows.
            if(!convertWindow(window, currentSampleConvertedWindows, 256, 30))
            {
                std::cerr << "[PopDel] Error in profile \"" << params.inputFiles[i] << "\"." << std::endl;
                std::ostringstream msg;
                msg << "[PopDel] Could not convert window " << params.contigNames[contigIdx][window.chrom] <<
                       ":" << window.beginPos << ". The profile \"" << params.inputFiles[i] << "\" might be corrupted.";
                SEQAN_THROW(IOError(toCString(msg.str())));
            }

            // Compare the read window with the current minimum.
            currentCoord.i1 = params.contigNames[contigIdx][currentSampleConvertedWindows[0].chrom];
            currentCoord.i2 = currentSampleConvertedWindows[0].beginPos;
            if (minCoord.i1 == "" || !lowerCoord(minCoord, currentCoord, params.nextRoi, params.allRois))
                minCoord = currentCoord;
        }
        // If there were not entries, try the next ROI.
        if (minCoord.i1 == params.nextRoi->seqName)
            break;
        else
            advanceRoiIterator(params.nextRoi, params.allRois); // TODO: Check potential for segfault
    }
    if (minCoord.i1 == "")
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not find any entries for any of the specified ROIs or reference sequence of" <<
               " the first sample. Terminating.";
        SEQAN_THROW(IOError(toCString(msg.str())));
    }
    std::ostringstream msg;
    msg << "The first window of all profiles starts at \'" << minCoord.i1 << ":" << minCoord.i2 << "\'.";
    printStatus(msg);
    return minCoord;
}
// ======================================================================================
// Function getFirstWindowInNextROI()
// ======================================================================================
// Return the first (=smallest) window position from all input streams in the next ROI.
// If there are no reads for the ROI in any sample, try the next ROI.
// Return false if there are no reads for any of the remaining ROI's, true otherwise.
inline bool getFirstWindowOnNextROI (Pair<CharString, __uint32> & minCoord,
                                     String<String<Window> >& convertedWindows,
                                     String<bool> & finishedROIs,
                                     PopDelCallParameters & params)
{
    minCoord = Pair<CharString, __uint32>("", maxValue<__uint32>());
    Pair<CharString, __uint32> currentCoord;
    std::string str;
    Window window;
    while (true)
    {
        for (unsigned i = 0; i < params.fileCount; ++i)
        {
            unsigned contigIdx = params.representativeContigs?0:i;
            String<Window> & currentSampleConvertedWindows = convertedWindows[i];
            std::ifstream in(toCString(params.inputFiles[i]), std::ios::in | std::ios::binary);
            if (!in.good())
                SEQAN_THROW(FileOpenError(toCString(params.inputFiles[i])));

            //jump to the region of interest.
            goNextRoi(in, i, params.nextRoi, finishedROIs, params);

            // decompress and read the window.
            readWindow(in, window, length(params.rgs[i]), params.uncompressedIn);
            if (isEmpty(window))
                continue;

            //Convert the 256bp window to a string of 30bp windows.
            if (!convertWindow(window, currentSampleConvertedWindows, 256, 30))
            {
                std::cerr << "[PopDel] Error in profile \"" << params.inputFiles[i] << "\"." << std::endl;
                std::ostringstream msg;
                msg << "[PopDel] Could not convert window at " << params.contigNames[contigIdx][window.chrom] <<
                ":" << window.beginPos << ". The profile \"" <<  params.inputFiles[i] << "\" might be corrupted.";
                SEQAN_THROW(IOError(toCString(msg.str())));
            }

            // Compare the read window with the current minimum.
            currentCoord.i1 = params.contigNames[contigIdx][currentSampleConvertedWindows[0].chrom];
            currentCoord.i2 = currentSampleConvertedWindows[0].beginPos;
            if (minCoord.i1 == "" || !lowerCoord(minCoord, currentCoord, params.nextRoi, params.allRois))
                minCoord = currentCoord;
        }
        // If there were not entries, try the next ROI.
        if (minCoord.i1 == params.nextRoi->seqName)
            break;
        else
        {
            if (peekNextRoiIterator(params.nextRoi, params.allRois))
                advanceRoiIterator(params.nextRoi, params.allRois);
            else
                return false;
        }
    }
    return true;
}
// =======================================================================================
// Function checkAndSwitch())
// =======================================================================================
// Check if the given beginPos is suited for the next subset of the chromosome profile and perform the
// neccessary switches otherwise.
// Return true, if the next segment can be loaded and false if loading has to be postponed.
inline bool checkAndSwitch(ChromosomeProfile & profile, const TReadGroupIndices & rg, const unsigned beginPos)
{
    SEQAN_ASSERT(!empty(rg));
    if (profile.startProfiles[rg[0]].tooBigForNext(beginPos))
    {
//         bool allEmpty = true;
//         for (unsigned i = 0; i < length(profile.startProfiles); ++i)
//         {
//             allEmpty &= empty(profile.startProfiles[i].sets[profile.startProfiles[i].writeSet].subset);
//             allEmpty &= empty(profile.endProfiles[i].sets[profile.endProfiles[i].writeSet].subset);
//             allEmpty &= profile.activeReads[i].empty();
//             if (!allEmpty)
//                 break;
//         }
//         if (allEmpty)
//         {
//             profile.resetTo(beginPos);   //TODO: We must NOT use beginPos, but the min pos off all next reads of all samples!
//             return true;                 // Re-implement this for significant speed-up
//         }
//         else
        {
            performSwitches(profile, rg);
            return false;
        }
    }
    else if(profile.startProfiles[rg[0]].needsSwitch(beginPos))           // Pos is too big for the current subset.
    {
        performSwitches(profile, rg);
        return false;
    }
    if (profile.endProfiles[rg[0]].readSet == profile.endProfiles[rg[0]].writeSet)
    {
        for (unsigned r = 0; r < length(rg); ++r)
        {
            profile.endProfiles[rg[r]].switchWriteSet();
            profile.endProfiles[rg[r]].correctConsecutiveSwitch(profile.activeLoad[rg[r]]);
        }
    }
    return true;
}
// =======================================================================================
// Function addRgRecordsToProfile()
// =======================================================================================
// Add all records of the read group to the ChromosomeProfile.
inline void addRgRecordsToProfile(ChromosomeProfile & profile,
                                  const TReadGroupIndices & rg,
                                  const String<Histogram> & histograms,
                                  const Window & win)
{
    for (unsigned r = 0; r < length(win.insertSizes); ++r)
    {   // Add all records of this read group to the profile.
        for (unsigned i = 0; i < length(win.insertSizes[r]); ++i)
        {
            const Pair<unsigned> & record = win.insertSizes[r][i];
            int innerDist = static_cast<int32_t>(record.i2) + histograms[rg[r]].median - 2 * histograms[rg[r]].readLength;
            unsigned endPos = record.i1;
            if (innerDist > 0)
                endPos += innerDist;
            profile.add(rg[r], record.i1, endPos, static_cast<int32_t>(record.i2));
        }
    }
}
// =======================================================================================
// Function readTillRoi()
// =======================================================================================
// Read all entries until the beginning of the ROI is reached, ignoring everything befor the ROI.
// Perform the necessary switches if necessary.
// Return 0 on success, 2 if the contig is not present (or lies before the starting point of reading) and 3 at EOF.
template<typename TStream>
inline unsigned readTillRoi(ChromosomeProfile & profile,
                            TStream & file,
                            Window & window,
                            const String<CharString> & contigNames,
                            Pair<CharString, __uint32> & coordinate,
                            const TReadGroupIndices & rg,
                            const GenomicRegion & roi)
{
    do
    {
        if (!readWindow(file, window, length(rg)))
        {
            performPartialSwitches(profile, rg);
            coordinate.i1 = "";
            coordinate.i2 = maxValue<unsigned>();
            return 3;               // EOF
        }
        if (contigNames[window.chrom] != coordinate.i1)
        {
            performPartialSwitches(profile, rg);
            coordinate.i1 = contigNames[window.chrom];
            coordinate.i2 = window.beginPos;
            return 2;               // Contig not present in sample or passed.
        }
    }
    while (window.beginPos + 255 < roi.beginPos);
    return 0;
}
// =======================================================================================
// Function readSegment()
// =======================================================================================
// Read all entries until the subset is full or the end of the contig or roi has been reached.
// Update coordinate to point to the first position at the next segment of the sample.
// Return 0 if the segment is full.
// Return 1 if the desired chromosome has not yet been reached.
// Return 2 if the end of the chromosome or ROI has been reached.
// Return 3 if EOF has been reached.
template<typename TStream>
inline unsigned readSegment(ChromosomeProfile & profile,
                            TStream & file,
                            const CharString & fileName,
                            const TReadGroupIndices & rg,
                            const String<Histogram> & histograms,
                            const String<CharString> & contigNames,
                            Pair<CharString, __uint32> & coordinate,
                            const GenomicRegion & roi,
                            String<Window> & convertedWindows)
{
    Window window;
    while (true)
    {
        unsigned ret = readTillRoi(profile, file, window, contigNames, coordinate, rg, roi);
        if (ret > 0)
            return ret;

        if(!convertWindow(window, convertedWindows, 256, 30))
        {
            std::cerr << "[PopDel] Error in profile \"" << fileName << "\"." << std::endl;
            std::ostringstream msg;
            msg << "[PopDel] Could not convert window at " << contigNames[window.chrom] <<
            ":" << window.beginPos << ". The profile \"" << fileName << "\" might be corrupted.";
            SEQAN_THROW(IOError(toCString(msg.str())));
        }
        for (Iterator<const String<Window> >::Type it = begin(convertedWindows); it != end(convertedWindows); ++it)
        {
            if (contigNames[it->chrom] != roi.seqName || it->beginPos >= roi.endPos)
            {
                performPartialSwitches(profile, rg);
                coordinate.i1 = contigNames[it->chrom];
                coordinate.i2 = it->beginPos;
                return 2;           // Contig has been completely processed.
            }
            if (it->beginPos < roi.beginPos)
                continue;

            if (profile.startProfiles[rg[0]].tooBigForNext(it->beginPos))
            {
                if (checkAllEmpty(profile))
                {
                    profile.resetTo(it->beginPos);
                }
                else
                {
                    performPartialSwitches(profile, rg);
                    coordinate.i1 = contigNames[it->chrom];
                    coordinate.i2 = it->beginPos;
                    return 0;
                }
            }
            else if(profile.startProfiles[rg[0]].needsSwitch(it->beginPos))           // Pos is too big for the current subset.
            {
                performPartialSwitches(profile, rg);
                coordinate.i1 = contigNames[it->chrom];
                coordinate.i2 = it->beginPos;
                return 0;
            }
            addRgRecordsToProfile(profile, rg, histograms, *it);
        }
    }
}
// Wrapper for handling de-compression if necessary
inline unsigned readSegment(ChromosomeProfile & profile,
                            std::ifstream & file,
                            const CharString & fileName,
                            const TReadGroupIndices & rg,
                            const String<Histogram> & histograms,
                            const String<CharString> & contigNames,
                            Pair<CharString, __uint32> & coordinate,
                            const GenomicRegion & roi,
                            String<Window> & convertedWindows,
                            const bool & uncompressed)
{
    if (uncompressed)
    {
        return readSegment(profile,
                           file,
                           fileName,
                           rg,
                           histograms,
                           contigNames,
                           coordinate,
                           roi,
                           convertedWindows);
    }
    else
    {
        zlib_stream::zip_istream unzipper(file);
        return readSegment(profile,
                           unzipper,
                           fileName,
                           rg,
                           histograms,
                           contigNames,
                           coordinate,
                           roi,
                           convertedWindows);
    }
}
// =======================================================================================
// Function processSegmentCode()
// =======================================================================================
// Process the code produced by readSegment()
// Set the nextReadPos and nextChromWindow accordingly.
// Also take care of setting the bools for finished files and ROIs.
// Return true if the processing of the file may proceed, false otherwise.
inline bool processSegmentCode(unsigned & nextReadPos,
                               String<bool> & finishedROIs,
                               String<bool> & finishedFiles,
                               PopDelCallParameters & params,
                               Windows & windows,
                               const unsigned & i,
                               const String<Pair<CharString, uint32_t> > & nextCandidateWindows,
                               const unsigned & segmentCode)
{
    if (segmentCode == 0u)
    {
        if (nextCandidateWindows[i].i2 < nextReadPos)
            nextReadPos = nextCandidateWindows[i].i2;
        return false;
    }
    else if (segmentCode < 3u)             // End of chromsome/ROI
    {
        finishedROIs[i] = true;
        if (windows.nextChromWindow.i1 == "" ||
            lowerCoord(nextCandidateWindows[i], windows.nextChromWindow, params.nextRoi, params.allRois))
        {
            windows.nextChromWindow = nextCandidateWindows[i];
        }
        return true;
    }
    else                                  //  segement code == 3 -> EOF
    {
        finishedROIs[i] = true;
        finishedFiles[i] = true;
        SEQAN_ASSERT_GT(params.fileCount, 0u);
        --params.fileCount;
        return true;
    }
}
#endif /* LOAD_PROFILE_POPDEL_CALL_H_ */
