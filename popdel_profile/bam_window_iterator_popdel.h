#ifndef BAM_WINDOW_ITERATOR_POPDEL_H_
#define BAM_WINDOW_ITERATOR_POPDEL_H_

#include "window_podel.h"
#include "bam_qual_req_popdel.h"
#include "profile_parameter_parsing_popdel.h"
#include "../utils_popdel.h"

using namespace seqan;

// =======================================================================================
// Struct BamWindowIterator
// =======================================================================================
// Actually not only for BAM files, but also for CRAM files.
struct BamWindowIterator
{
    typedef std::map<CharString, unsigned> TReadGroups;    // Mapping of ReadGroupID:Index
    HtsFileIn infile;                                      // The input BAM/CRAM-file.
    TReadGroups readGroups;                                // Map of read group IDs and their rank (order of apperance).
    BamQualReq qualReqFwd;                                 // Quality requirements for forward reads.
    BamQualReq qualReqRev;                                 // Quality requirements for reverse reads.
    std::map<CharString, Pair<__int32> > goodFwdReads;     // Forward reads that meet all quality requirements;
                                                               // stored as: map<ID, <innerBeginPos, insertBeginPos> >.
    BamAlignmentRecord nextRecord;                         // The next record to be processed.
    __int32 nextRecordsInnerBegin;                         // Begin of the next records's inner part (part betw. reads).
    __int32 nextRecordsInsertBegin;                        // Begin of next record's insert.
    __int32 maxDeletionSize;                               // Biggest insertSize of all RGs to be considered.
    String<int32_t> maxInsertSizes;                        // Max. insertSizes per RG.
    __int32 windowSize;                                    // Size of one window.
    __int32 windowShift;                                   // Window shift per iteration (=window size by default).
    String<Window> activeWindows;                          // String of currently active windows.
    __int32 currentWindowBegin;                            // Begin position of current Window.
    String<GenomicRegion> intervals;                       // Intervals to work on.
    Iterator<String<GenomicRegion> >::Type currentReadingInterval; //
    Iterator<String<GenomicRegion> >::Type currentWindowInterval;  //
    bool mergeRG;

    BamWindowIterator(const CharString & filename, const CharString & referenceFile)
    {
        if(!open(infile, toCString(filename), referenceFile == "" ? NULL : toCString(referenceFile)))
            SEQAN_THROW(IOError("Could not open sequence file."));

        if (infile.fp->fn_aux)
        {
            if (infile.fp->is_cram)
            {
                std::ostringstream msg;
                msg << "Using \'" << referenceFile << "\' as reference file instead of header entries.";
                printStatus(msg);
            }
            else
            {
                std::ostringstream msg;
                msg << "Warning: Sequence file is not in CRAM format but a reference file has been specified. The reference file will be ignored.";
                printStatus(msg);
            }
        }

        if(!loadBaiCrai(infile))
            SEQAN_THROW(IOError("Could not open index file."));
    }
};
// ---------------------------------------------------------------------------------------
// Function posToWindow()
// ---------------------------------------------------------------------------------------
// Returns an iterator to the window which overlaps the position 'pos'.
inline Iterator<String<Window> >::Type posToWindow(const unsigned & pos, BamWindowIterator & bwi)
{
    unsigned windowIndex =  pos / bwi.windowShift;          // Begin position of the window that overlaps pos.
    unsigned w = windowIndex % length(bwi.activeWindows);   // Index of this window in the String of windows.
    return begin(bwi.activeWindows) + w;                    // Return iterator to w'th window.
}
// =======================================================================================
// Function operator*()
// =======================================================================================
// Defines the de-reference operator for the BamWindowIterator, thus returning the active window.
inline Reference<String<Window> >::Type operator*(BamWindowIterator & bwi)
{
    return *posToWindow(bwi.currentWindowBegin, bwi);
}
// ---------------------------------------------------------------------------------------
// Function processRecord()
// ---------------------------------------------------------------------------------------
// Return 1 if processing of read (nextRecord) is done (=Checking requirements and maybe adding its info to bwi).
// Return 0 if a good rev. read with a good mate was found whose processing needs to wait until window is reached.
// Return -1 if the record lies outside of the interval specified in bwi.interval.
inline int processRecord(BamWindowIterator & bwi)
{
    // Process forward reads.
    if (!isReverseRead(bwi.nextRecord))                            // Check if read is first in pair.
    {
        if (bwi.nextRecord.rID == bwi.currentReadingInterval->rID &&
            bwi.nextRecord.beginPos >= bwi.currentReadingInterval->endPos)
        {
            return -1;                                                      // Fwd.-read starts after current interval.
        }
        if (bwi.nextRecord.pNext < bwi.currentReadingInterval->beginPos)    // Rev.-read ends before current interval.
        {
            return 1;
        }
        if (meetsRequirements(bwi.nextRecord, bwi.qualReqFwd) && bwi.nextRecord.beginPos <= bwi.nextRecord.pNext)
        {
            __int32 innerBeginPos = bwi.nextRecord.beginPos + getAlignmentLengthInRef(bwi.nextRecord);
            __int32 insertBeginPos = bwi.nextRecord.beginPos;
            if (bwi.nextRecord.cigar[0].operation == 'S')                   // If next record is soft-clipped...
            {
                insertBeginPos -= bwi.nextRecord.cigar[0].count;            // ...adjust its starting position.
            }
            if (bwi.nextRecord.tLen <= bwi.maxInsertSizes[getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG)])
            {
                bwi.goodFwdReads[bwi.nextRecord.qName] = Pair<__int32>(innerBeginPos, insertBeginPos);
            }
            return 1;
        }
    }
    // Process reverse reads.
    else if (bwi.nextRecord.beginPos >= bwi.nextRecord.pNext)       // Check in case there is an inversion.
    {
        std::map<CharString, Pair<__int32> >::iterator fwd = bwi.goodFwdReads.find(bwi.nextRecord.qName);
        if (fwd != bwi.goodFwdReads.end()) // If the corresponding forward read is in the list of good forward reads...
        {   // ...save its innerBeginPos and insertBeginPos for later use in nextRecord and erase it from goodFwdReads.
            bwi.nextRecordsInnerBegin = (fwd->second).i1;
            bwi.nextRecordsInsertBegin = (fwd->second).i2;
            bwi.goodFwdReads.erase(fwd);
            if (meetsRequirements(bwi.nextRecord, bwi.qualReqRev))
                return 0;
        }
    }
    return true;
}
// ---------------------------------------------------------------------------------------
// Function tryReadRecord()
// ---------------------------------------------------------------------------------------
// Try to read the next record in the set region. Do not read the sequence or qualities.
// Return true on sucess and false if EOF has been reached.
inline bool tryReadRecord(BamAlignmentRecord & record, HtsFileIn & infile)
{
            bool ret = false;
            try
            {
                ret = seqFreeReadRegion(record, infile);
            }
            catch (Exception const & e)
            {
                SEQAN_THROW(IOError("Could no read record in BAM-File."));
            }
            return ret;
}
// ---------------------------------------------------------------------------------------
// Function goToInterval()
// ---------------------------------------------------------------------------------------
// Uses the Bam index to jump to the region of interest and reads the records until the specified interval is reached.
inline void goToInterval(BamWindowIterator & bwi)
{
    GenomicRegion const & itv = *bwi.currentReadingInterval;
    setRegion(bwi.infile, toCString(itv.seqName), itv.beginPos + 4, itv.endPos);
    bool notAtEnd = true;
    do
    {
        notAtEnd = tryReadRecord(bwi.nextRecord, bwi.infile);
    }
    while (bwi.nextRecord.rID == itv.rID && bwi.nextRecord.beginPos < itv.beginPos && notAtEnd);
}
// ---------------------------------------------------------------------------------------
// Function goToNextReverseRecord()
// ---------------------------------------------------------------------------------------
// Find the next reverse record in the intervals and process it using _processRecord. Called by _goToFirstReverseRecord.
// Return true if a good rev. read with a good mate was found whose processing needs to wait until window is reached.
// Return false if end of file is reached.
inline bool goToNextReverseRecord(BamWindowIterator & bwi)
{

    while(bwi.currentReadingInterval != end(bwi.intervals))
    {
        while (tryReadRecord(bwi.nextRecord, bwi.infile))
        {

            int pr = processRecord(bwi);
            if (pr == 0)                                        // Only good reverse records cannot be completely processed.
            {
                return true;
            }
            else if (pr < 0 && bwi.goodFwdReads.empty())        // Record is outside of interval.
            {   // TODO Check if this can still happen.
                ++bwi.currentReadingInterval;
                if (bwi.currentReadingInterval >= end(bwi.intervals))
                {
                    return false;                               // Last interval has been processed.
                }
                goToInterval(bwi);
            }
        }
        if (bwi.currentReadingInterval != end(bwi.intervals)) // Check the next interval
        {
            ++bwi.currentReadingInterval;
            if (bwi.currentReadingInterval >= end(bwi.intervals))
            {
                return false;                               // Last interval has been processed.
            }
            goToInterval(bwi);
        }
    }
    return false;                                           // End of file reached.
}
// ---------------------------------------------------------------------------------------
// Function goToFirstReverseRecord()
// ---------------------------------------------------------------------------------------
// Performs the same task as _goToNextReverseRecord, except that it also takes care of the intizializations.
// Should be called for the first record in the first interval.
// Return true if a good rev. read with a good mate was found whose processing needs to wait until window is reached.
// Return false end of file is reached.
// Call _goToNextReverseRecord, if no good reverse record can be found in the first try.
inline bool goToFirstReverseRecord(BamWindowIterator & bwi)
{
    setRIDs(bwi.intervals, bwi.infile);
    // Initialize the first interval and read the first record after the interval's begin from infile
    bwi.currentReadingInterval = begin(bwi.intervals);
    goToInterval(bwi);
    while (bwi.nextRecord.rID > (*bwi.currentReadingInterval).rID)  // Is record in current interval?
    {
        ++bwi.currentReadingInterval;
        if (bwi.currentReadingInterval >= end(bwi.intervals))
        {
            return false;                                   // Last interval has been processed.
        }
        goToInterval(bwi);
    }
    int pr = processRecord(bwi);
    if (pr == 0)                                            // Only good reverse records cannot be completely processed.
    {
        return true;
    }
    else if (pr < 0 && bwi.goodFwdReads.empty())            // Record is outside of interval.
    {
        ++bwi.currentReadingInterval;
        if (bwi.currentReadingInterval >= end(bwi.intervals))
        {
            return false;                               // Last interval has been processed.
        }
        goToInterval(bwi);
    }
    return goToNextReverseRecord(bwi);     // Continue searching via _goToNextReverseRecord.
}
// ---------------------------------------------------------------------------------------
// Function addRecordToWindow()
// ---------------------------------------------------------------------------------------
// Add information about the number of overlapped windows and the insert size of the record to the window of read group.
inline void addRecordToWindows(BamWindowIterator & bwi)
{
    typedef Iterator<String<Window> >::Type TWindowIter;
    TWindowIter it = posToWindow(bwi.nextRecordsInnerBegin, bwi);
    TWindowIter itEnd = posToWindow(bwi.nextRecord.beginPos, bwi);
    TWindowIter wEnd = end(bwi.activeWindows);
    itEnd++;
    if (itEnd == wEnd)
        itEnd = begin(bwi.activeWindows);
    unsigned insertSize = bwi.nextRecord.beginPos + getAlignmentLengthInRef(bwi.nextRecord) - bwi.nextRecordsInsertBegin;
    if (bwi.nextRecord.cigar[length(bwi.nextRecord.cigar) - 1].operation == 'S')    // If next record is soft-clipped...
        insertSize += bwi.nextRecord.cigar[length(bwi.nextRecord.cigar) - 1].count; //...adjust insert size
    /*unsigned numWindows = 1u;                                                  // Number of windows the record overlaps.
    if (bwi.nextRecordsInnerBegin < bwi.nextRecord.beginPos)
    {
        if (it < itEnd)
            numWindows = itEnd - it;
        else
            numWindows = wEnd - it + itEnd - begin(bwi.activeWindows);
    }*/
    unsigned rg = getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG);
    addRecord(*it, bwi.nextRecordsInnerBegin, insertSize, rg, length(bwi.readGroups));
}
// ---------------------------------------------------------------------------------------
// Function resetWindows()
// ---------------------------------------------------------------------------------------
// Reset all active windows in the current interval.
inline void resetWindows(BamWindowIterator & bwi)
{
    SEQAN_ASSERT_LT(bwi.currentReadingInterval, end(bwi.intervals));
    // Set the interval to the next interval with good read pairs.
    bwi.currentWindowInterval = bwi.currentReadingInterval;
    // Set the window begin position.
    int maxInsertSize = max(bwi.maxInsertSizes);
    //int maxInsertSize = bwi.maxInsertSizes[getReadGroup(bwi.nextRecord.tags, bwi.readGroups)];
    unsigned pos = bwi.nextRecord.beginPos > maxInsertSize ? bwi.nextRecord.beginPos - maxInsertSize : 0;
    bwi.currentWindowBegin = (pos / bwi.windowShift) * bwi.windowShift;
    // Clear all windows and set their chromosome and begin position.
    for (__int32 i = bwi.currentWindowBegin;
         i < bwi.currentWindowBegin + maxInsertSize + bwi.windowSize + 1;
         i += bwi.windowShift)
        *posToWindow(i, bwi) = Window((*bwi.currentWindowInterval).rID, i);
}
// ---------------------------------------------------------------------------------------
// Function isLastWindow()
// ---------------------------------------------------------------------------------------
// Return true if the last window of the (globally) last interval has been processed, false otherwise.
inline bool isLastWindow(BamWindowIterator & bwi)
{
    if (bwi.currentReadingInterval >= end(bwi.intervals) &&
        bwi.currentWindowBegin >= (*bwi.currentWindowInterval).endPos)
    {
        std::ostringstream msg;
        msg << "Finished scanning interval \'" << (*bwi.currentWindowInterval).seqName << ":";
        msg << (*bwi.currentWindowInterval).beginPos + 1 << "-" << (*bwi.currentWindowInterval).endPos << "\'.";
        printStatus(msg);
        return true;
    }
    return false;
}
// ---------------------------------------------------------------------------------------
// Function isLastWindowOnSeq()
// ---------------------------------------------------------------------------------------
//Return true if the last window of an interval/sequence/contig/chromosome has been processed, false otherwise.
inline bool isLastWindowOnSeq(BamWindowIterator & bwi)
{
    SEQAN_ASSERT(!isLastWindow(bwi));
    if (bwi.currentWindowBegin > (*bwi.currentWindowInterval).endPos)
    {
        std::ostringstream msg;
        msg << "Finished scanning interval \'" << (*bwi.currentWindowInterval).seqName << ":";
        msg << (*bwi.currentWindowInterval).beginPos + 1 << "-" << (*bwi.currentWindowInterval).endPos << "\'.";
        printStatus(msg);
        return true;
    }
    return false;
}
// =======================================================================================
// Function goNext()
// =======================================================================================
// Shift the iterator to the next window. Return true on success, false if the last window has been reached.
inline bool goNext(BamWindowIterator & bwi)
{
    while (true)
    {
        // Reset the previously reported/skipped window.
        *bwi = Window((*bwi.currentWindowInterval).rID,
                      bwi.currentWindowBegin + bwi.windowShift * length(bwi.activeWindows));
        // Move active windows by one.
        bwi.currentWindowBegin += bwi.windowShift;
        // Return if we reached last window in last interval.
        if (isLastWindow(bwi))
            return false;
        // Reset all active windows if we reached last window in an interval.
        if (isLastWindowOnSeq(bwi))
            resetWindows(bwi);
        // Add all records that may overlap the first active window.
        int maxInsertSize = max(bwi.maxInsertSizes);
        //int maxInsertSize = bwi.maxInsertSizes[getReadGroup(bwi.nextRecord.tags, bwi.readGroups)];
        while (bwi.currentReadingInterval == bwi.currentWindowInterval &&
               bwi.nextRecordsInnerBegin - maxInsertSize < bwi.currentWindowBegin + bwi.windowSize)
        {
            addRecordToWindows(bwi);
            goToNextReverseRecord(bwi);
        }
        // Return if the first active window is not empty.
        if (!isEmpty(*bwi))
            return true;
    }
}
// =======================================================================================
// Function calculateMaxInsertSizes()
// =======================================================================================
// Sets the calculateMaxInsertSizes for each read group.
// Return the biggest maxInsertSize of all read groups.
inline void calculateMaxInsertSizes(BamWindowIterator & bwi,
                                    const String<Histogram> & hists)
{
    unsigned rgCount = length(hists);
    resize(bwi.maxInsertSizes, rgCount);
    for (unsigned rg = 0; rg < rgCount; ++rg)
        bwi.maxInsertSizes[rg] = hists[rg].median + 2 * hists[rg].stddev + bwi.maxDeletionSize; // TODO: Check!
}
// =======================================================================================
// Function initialize()
// =======================================================================================
// Initialize all parameters for the window iterator, open BAM file and set window to the first good interval.
void initialize(BamWindowIterator & bwi,
                std::map<CharString, unsigned> & readGroups,
                String<GenomicRegion> & intervals,
                BamQualReq & qualReq,
                int maxDeletionSize,
                int windowSize,
                int windowShift,
                const String<Histogram> & hists,
                const bool mergeRG)
{
    // Set the parameters.
    bwi.mergeRG = mergeRG;
    bwi.windowSize = windowSize;
    bwi.windowShift = windowShift;
    if (mergeRG)
    {
        bwi.readGroups[readGroups.begin()->first] = readGroups.begin()->second;
    }
    else
    {
        bwi.readGroups = readGroups;
    }
    bwi.maxDeletionSize = maxDeletionSize;
    calculateMaxInsertSizes(bwi, hists);
    bwi.intervals = intervals;
    // Set the quality requirements for reads.
    bwi.qualReqFwd = qualReq;
    bwi.qualReqRev = qualReq;
    bwi.qualReqFwd.flagsSet |= BAM_FLAG_NEXT_RC;
    bwi.qualReqFwd.flagsUnset |= BAM_FLAG_RC;
    bwi.qualReqFwd.flagsSet &= ~BAM_FLAG_RC;
    bwi.qualReqFwd.flagsUnset &= ~BAM_FLAG_NEXT_RC;
    bwi.qualReqRev.flagsSet |= BAM_FLAG_RC;
    bwi.qualReqRev.flagsUnset |= BAM_FLAG_NEXT_RC;
    bwi.qualReqRev.flagsSet &= ~BAM_FLAG_NEXT_RC;
    bwi.qualReqRev.flagsUnset &= ~BAM_FLAG_RC;
    // Open the input bam file, read the header and go to the reverse read of the first good read pair in infile.
    if (!goToFirstReverseRecord(bwi))
        SEQAN_THROW(IOError("[PopDel] No BAM records meet the requirements."));
    // Initialize windows.
    int maxInsertSize = max(bwi.maxInsertSizes);    //TODO: Check if this works out and replace max with parameter.
    //int maxInsertSize = bwi.maxInsertSizes[getReadGroup(bwi.nextRecord.tags, bwi.readGroups)];
    resize(bwi.activeWindows, ((maxInsertSize + bwi.windowSize) / bwi.windowShift) + 1);
    resetWindows(bwi);
    // Add all records that may overlap the first window.
    while (bwi.currentReadingInterval == bwi.currentWindowInterval &&
        bwi.nextRecordsInnerBegin - maxInsertSize < bwi.currentWindowBegin + bwi.windowSize)
    {
        addRecordToWindows(bwi);
        goToNextReverseRecord(bwi);
    }
    // Go to next window if first window is empty.
    if (isEmpty(*bwi))
        goNext(bwi);
    std::ostringstream msg;
    msg << "Initialized window iterator with windowSize=" << bwi.windowSize << " and windowShift=" << windowShift << ".";
    printStatus(msg);
}

#endif // BAM_WINDOW_ITERATOR_POPDEL_H_
