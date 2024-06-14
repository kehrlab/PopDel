#ifndef BAM_WINDOW_ITERATOR_POPDEL_H_
#define BAM_WINDOW_ITERATOR_POPDEL_H_

#include "window_popdel.h"
#include "bam_qual_req_popdel.h"
#include "profile_parameter_parsing_popdel.h"
#include "profile_translocation_popdel.h"
#include "profile_translocation_popdel.h"
#include "../utils_popdel.h"

using namespace seqan;


struct GoodReadBufferEntry; // Forward declare for use in BamWindowInterator.

// =======================================================================================
// Struct BamWindowIterator
// =======================================================================================
// Actually not only for BAM files, but also for CRAM files.
struct BamWindowIterator
{
    typedef std::map<CharString, unsigned> TReadGroups;    // Mapping of ReadGroupID:Index
    HtsFileIn infile;                                      // The input BAM/CRAM-file.
    TReadGroups readGroups;                                // Map of read group IDs and their rank (order of apperance).
    BamQualReq qualReq;                                    // Quality requirements for forward reads.
    std::map<CharString, GoodReadBufferEntry> goodLeftReads;    // Left reads that meet all quality requirements;
    String<TranslocationBuffer> translocationBuffers;     // Buffers holding all pairs affected by translocations. 1 per RG.
    BamAlignmentRecord nextRecord;                        // The next record to be processed.
    int32_t leftTipPos;                                   // clipped 3' end of the leftmost read of the pair
    int32_t distance;                                     // FR-reads: 5'-end distance. Other reads: 3'-end distance.
    int32_t totalClipping;                                // Total ammount of soft clipping in both reads of the pair.
    std::vector<uint8_t> clipping;                        // Number of clipped bases at all read ends: (5', 3', 3', 5')
    Orientation orientation;                              // Orientation of the current pair
    int32_t maxDeletionSize;                              // Biggest tip distance of all RGs to be considered.
    String<int32_t> maxTipDistances;                      // Max. tip distance per RG.
    String<int32_t> readLengths;                          // Read length for every RG
    int32_t windowSize;                                   // Size of one window.
    int32_t windowShift;                                  // Window shift per iteration (=window size by default).
    String<Window> activeWindows;                         // String of currently active windows.
    int32_t currentWindowBegin;                           // Begin position of current Window.
    String<GenomicRegion> intervals;                      // Intervals to work on.
    TranslocationBlacklist blacklist;                      // Blacklist for translocation reads
    Iterator<String<GenomicRegion> >::Type currentReadingInterval;
    Iterator<String<GenomicRegion> >::Type currentWindowInterval;
    String<bool> finishedRIDs;                         // Bool for each rID in the BAM file indicating finished status.
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
// =======================================================================================
// Function initContigStatusList()
// =======================================================================================
// Prepare BamWindowIterator::finishedRIDs for later use
inline void initContigStatusList(BamWindowIterator & bwi, const unsigned numContigs)
{
    resize(bwi.finishedRIDs, numContigs, false, Exact());
}
// =======================================================================================
// Function getNextReadLenght()
// =======================================================================================
// Return the read lenght of the read group bwi.nextRead belongs to.
inline int getNextReadLenght(const BamWindowIterator & bwi)
{
    return bwi.readLengths[getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG)];
}
// =======================================================================================
// Function getLeftClip()
// =======================================================================================
// Return the number of sofclipped bases at the left end of the record.
inline uint32_t getLeftClip(const BamAlignmentRecord & rec)
{
    SEQAN_ASSERT_GT(length(rec.cigar), 0u);
    if (rec.cigar[0].operation == 'S')
        return rec.cigar[0].count;
    else
        return 0;
}
// =======================================================================================
// Function getLeftClip()
// =======================================================================================
// Return the number of sofclipped bases at the left end of BamWindowInterator::nextRecord.
inline uint32_t getLeftClip(const BamWindowIterator & bwi)
{
    return getLeftClip(bwi.nextRecord);
}
// =======================================================================================
// Function getRightClip()
// =======================================================================================
// Return the number of sofclipped bases at the right end of the record.
inline uint32_t getRightClip(const BamAlignmentRecord & rec)
{
    SEQAN_ASSERT_GT(length(rec.cigar), 0u);
    unsigned i = length(rec.cigar) - 1;
    if (rec.cigar[i].operation == 'S')
        return rec.cigar[i].count;
    else
        return 0;
}
// Return the number of sofclipped bases at the right end of BamWindowInterator::nextRecord.
inline uint32_t getRightClip(const BamWindowIterator & bwi)
{
    return getRightClip(bwi.nextRecord);
}
// =======================================================================================
// Struct GoodReadBufferEntry
// =======================================================================================
struct GoodReadBufferEntry
{
    int32_t clippedLeftEnd;         // = beginPos
    int32_t fullLeftEnd;            // = beginPos - leftClip
    int32_t clippedRightEnd;        // = beginPos + getAlignmentLengthInRef - 1
    int32_t fullRightEnd;           // = clippedRightEnd + rightClip
    int32_t flag;

    GoodReadBufferEntry(): fullLeftEnd(0), clippedRightEnd(0), fullRightEnd(0), flag(0){};

    GoodReadBufferEntry(const int32_t fle, const int32_t cre, const int32_t fre, const int32_t flg):
    fullLeftEnd(fle),
    clippedRightEnd(cre),
    fullRightEnd(fre),
    flag(flg){};

    GoodReadBufferEntry(const BamWindowIterator & bwi)
    {
        clippedLeftEnd = bwi.nextRecord.beginPos;
        fullLeftEnd = clippedLeftEnd - getLeftClip(bwi);
        if (fullLeftEnd < 0)
            fullLeftEnd = 0;
        clippedRightEnd =  clippedLeftEnd + getAlignmentLengthInRef(bwi.nextRecord) - 1;
        fullRightEnd =  clippedRightEnd + getRightClip(bwi);
        flag = bwi.nextRecord.flag;
    };
};
// Overload for directly applying getPairOrientation() on a GoodReadBufferEntry.
inline Orientation getPairOrientation(const GoodReadBufferEntry & entry)
{
    return getPairOrientation(entry.flag);
}
// ---------------------------------------------------------------------------------------
// Function posToWindow()
// ---------------------------------------------------------------------------------------
// Return an iterator to the window which overlaps the position 'pos'.
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
// =======================================================================================
// Function isLeftAlignment()
// =======================================================================================
// Return true if the startin position of the given record is smaller than that of the mate
// Does not correct for softclipping!
inline bool isLeftAlignment(const BamAlignmentRecord & record)
{
    return record.beginPos < record.pNext;
}
// =======================================================================================
// Function isRightAlignment()
// =======================================================================================
// Return true if the startin position of the given record is bigger than that of the mate
// Does not correct for softclipping!
inline bool isRightAlignment(const BamAlignmentRecord & record)
{
    return record.beginPos > record.pNext;
}
// =======================================================================================
// Function getAbsTipDistance()
// =======================================================================================
// Return the maximum distance of the clipped 3'ends to each other, i.e. the tip to tip distance.
// Because we tLen is ambigous, we cannot give the exact end of the second's read 3' end, without looking at it.
// Therefore, we calculate the highest possible distance between the tips, based on the read length and pNExt.
// For normally oriented read pairs this equals to the inner distance.
inline int32_t getAbsTipDistance(const BamWindowIterator & bwi)
{
    int32_t a = bwi.nextRecord.beginPos;
    if (!hasFlagRC(bwi.nextRecord))
        a += getAlignmentLengthInRef(bwi.nextRecord) - 1;

    int32_t b;
    if (hasFlagNextRC(bwi.nextRecord))
        b = bwi.nextRecord.pNext;
    else
        b = bwi.nextRecord.pNext + getNextReadLenght(bwi);
    return abs(a - b);
}
// =======================================================================================
// Function inMaxRange()
// =======================================================================================
// Return true if the 3' ends of the record are close enough to be added to the window.
// Return false otherwise.
inline bool inMaxRange(const BamWindowIterator & bwi)
{
    int rg = getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG);
    return getAbsTipDistance(bwi) + bwi.readLengths[rg] <= bwi.maxTipDistances[rg];
    // We need the addition of the read length to make sure that we don't allow read pairs that might
    // start after a window has already been closed, because we found another read with a lower left tip.
}
// =======================================================================================
// Function inWindowRange()
// =======================================================================================
// Return true if the left tip of bwi does lie in the valid range of windows.
// Return false otherwise.
inline bool inWindowRange(const BamWindowIterator & bwi, const int & maxTipDistance)
{
    return bwi.leftTipPos < bwi.currentWindowBegin + bwi.windowSize + maxTipDistance;
}
//
// =======================================================================================
// Function processLeftAlignment()
// =======================================================================================
// Processes a (probable) left alignment by adding its properties to the map of good reads.
// Return false if the read lies after the current interval (i.e the interval has been completely processed).
// Return true otherwise.
inline bool processLeftAlignment(BamWindowIterator & bwi)
{
    SEQAN_ASSERT(!isTranslocated(bwi.nextRecord));
    if (bwi.nextRecord.beginPos >= bwi.currentReadingInterval->endPos)
        return false;                                                   // Fwd. read starts after current interval.
    if (bwi.nextRecord.pNext < bwi.currentReadingInterval->beginPos)
        return true;                                                    // Rev. read ends before current interval.
    // if (inMaxRange(bwi))                                             // general-purpose profiles do not need max distance
    bwi.goodLeftReads[bwi.nextRecord.qName] = GoodReadBufferEntry(bwi);
    return true;
}
// =======================================================================================
// Function getTip()
// =======================================================================================
// Return the position of the record's 3'-end.
// This equals r.beginPos for reverse complement reads and r.beginPos + getAlignmentLengthInRef(r) -1 for forward reads.
inline int32_t getTip(const BamAlignmentRecord & r)
{
    if (hasFlagRC(r))
        return r.beginPos;
    else
        return r.beginPos + getAlignmentLengthInRef(r) - 1;
}
inline int32_t getTip(const GoodReadBufferEntry & r)
{
    if (r.flag & 16)         // read is reverse complement
        return r.clippedLeftEnd;
    else
        return r.clippedRightEnd;
}
// =======================================================================================
// Function assignTipAndDistance()
// =======================================================================================
// Assign the position of the tip of the left read and the distance between the tips (3' ends) of the reads.
inline void assignTipAndDistance(int32_t & leftTipPos,
                                 int32_t & distance,
                                 const GoodReadBufferEntry & leftRead,
                                 const BamAlignmentRecord & rightRead)
{
    leftTipPos = getTip(leftRead);
    distance = getTip(rightRead) - leftTipPos;
}
//Overload in case the first read is the BamAlignmentRecord record.
inline void assignTipAndDistance(int32_t & leftTipPos,
                                 int32_t & distance,
                                 const BamAlignmentRecord & leftRead,
                                 const GoodReadBufferEntry & rightRead)
{
    leftTipPos = getTip(leftRead);
    distance = getTip(rightRead) - leftTipPos;
}
// =======================================================================================
// Function getFullClip()
// =======================================================================================
// Return the sum of leftClip and rightClip of a read.
inline int getFullClip(const GoodReadBufferEntry & r)
{
    return (r.clippedLeftEnd - r.fullLeftEnd) + (r.fullRightEnd - r.clippedRightEnd);
}
inline unsigned getFullClip(const BamAlignmentRecord & r)
{
    return getLeftClip(r) + getRightClip(r);
}
// =======================================================================================
// Function getTotalClipping()
// =======================================================================================
// Return the sum of fullClip of both reads of a pair.
inline int getTotalClipping(const GoodReadBufferEntry & g,
                            const BamAlignmentRecord & b)
{
    return getFullClip(g) + getFullClip(b);
}
inline int getTotalClipping(const BamAlignmentRecord & b,
                            const GoodReadBufferEntry & g)
{
    return getTotalClipping(g, b);
}
// =======================================================================================
// Function getClipping()
// =======================================================================================
// Return the individual number of clipped bases at all read ends in a pair
// assumes that r1 is the first read in pair (according to alignment position)
inline std::vector<uint8_t> getClipping(const BamAlignmentRecord & r1, const BamAlignmentRecord & r2)
{
    std::vector<uint8_t> clipping{0, 0, 0, 0};
    if (hasFlagRC(r1)) // read is reverse complement
    {
        clipping[0] = getRightClip(r1);
        clipping[1] = getLeftClip(r1);
    } else {
        clipping[0] = getLeftClip(r1);
        clipping[1] = getRightClip(r1);
    }

    if (hasFlagRC(r2)) // read is reverse complement
    {
        clipping[2] = getLeftClip(r2);
        clipping[3] = getRightClip(r2);
    } else {
        clipping[2] = getRightClip(r2);
        clipping[3] = getLeftClip(r2);
    }
    return clipping;
}

inline std::vector<uint8_t> getClipping(const BamAlignmentRecord & r, const GoodReadBufferEntry & g)
{
    std::vector<uint8_t> clipping{0, 0, 0, 0};
    if (hasFlagRC(r)) // read is reverse complement
    {
        clipping[0] = getRightClip(r);
        clipping[1] = getLeftClip(r);
    } else {
        clipping[0] = getLeftClip(r);
        clipping[1] = getRightClip(r);
    }

    if (g.flag & 16) // read is reverse complement
    {
        clipping[2] = g.clippedLeftEnd - g.fullLeftEnd;
        clipping[3] = g.fullRightEnd - g.clippedRightEnd;
    } else {
        clipping[2] = g.fullRightEnd - g.clippedRightEnd;
        clipping[3] = g.clippedLeftEnd - g.fullLeftEnd;
    }
    return clipping;
}

inline std::vector<uint8_t> getClipping(const GoodReadBufferEntry & g, const BamAlignmentRecord & r)
{
    std::vector<uint8_t> clipping{0, 0, 0, 0};

    if (g.flag & 16) // read is reverse complement
    {
        clipping[0] = g.fullRightEnd - g.clippedRightEnd;
        clipping[1] = g.clippedLeftEnd - g.fullLeftEnd;
    } else {
        clipping[0] = g.clippedLeftEnd - g.fullLeftEnd;
        clipping[1] = g.fullRightEnd - g.clippedRightEnd;
    }

    if (hasFlagRC(r)) // read is reverse complement
    {
        clipping[2] = getLeftClip(r);
        clipping[3] = getRightClip(r);
    } else {
        clipping[2] = getRightClip(r);
        clipping[3] = getLeftClip(r);
    }
    return clipping;
}

// =======================================================================================
// Function processRightAlignment()
// =======================================================================================
// Processes a (probable) right alignment by correcting for the clipping and adding it to the map of good reads.
// Also checks if the left an right alignment have been confused due to sofclipping.
// Return true if a good right alignment has been processed can can be written to ouput, false otherwise.
inline bool processRightAlignment(BamWindowIterator & bwi)
{
    SEQAN_ASSERT(!isTranslocated(bwi.nextRecord));
    std::map<CharString, GoodReadBufferEntry>::iterator fwd = bwi.goodLeftReads.find(bwi.nextRecord.qName);
    if (fwd != bwi.goodLeftReads.end())
    {
        int32_t fullLeftEnd = bwi.nextRecord.beginPos - getLeftClip(bwi);
        if (fullLeftEnd < fwd->second.fullLeftEnd)
        {   // left and right read have been confused due to soft clipping. 'bwi.nextRecord' is the left read.
            bwi.orientation = getPairOrientation(fwd->second);
            assignTipAndDistance(bwi.leftTipPos, bwi.distance, bwi.nextRecord, fwd->second);
            bwi.clipping = getClipping(bwi.nextRecord, fwd->second);
        }
        else
        {
            bwi.orientation = getPairOrientation(bwi.nextRecord);
            assignTipAndDistance(bwi.leftTipPos, bwi.distance, fwd->second, bwi.nextRecord);
            bwi.clipping = getClipping(fwd->second, bwi.nextRecord);
        }
        bwi.totalClipping = getTotalClipping(bwi.nextRecord, fwd->second);
        SEQAN_ASSERT_GEQ(bwi.totalClipping, 0);
        bwi.goodLeftReads.erase(fwd);
        return true;
    }
    return false;
}
// =======================================================================================
// Function checkAndRemoveTranslocatedRecord()
// =======================================================================================
// Check the if the mate of a bad translocated record possibly is in the translocationBuffer
// and try to remove it if so.
inline void checkAndRemoveTranslocatedRecord(BamWindowIterator & bwi)
{
    if (rIDAlreadyProcessed(bwi.finishedRIDs, bwi.nextRecord.rNextId))
    {
        if (!bwi.blacklist.contains(bwi.nextRecord.rNextId, bwi.nextRecord.pNext))
        {
            unsigned rg = getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG);
            bwi.translocationBuffers[rg].removeIfFound(bwi.nextRecord.qName);
        }
    }
}
// =======================================================================================
// Function checkAndAddTranslocatedRecord()
// =======================================================================================
// Check if the mate of the record has already been added to the translocationBuffer
// and insert a new entry / update the existing one accordingly if both reads lie in ROIS and are not blacklisted.
inline void checkAndAddTranslocatedRecord(BamWindowIterator & bwi)
{
    if (rIDAlreadyProcessed(bwi.finishedRIDs, bwi.nextRecord.rNextId))
    {
        if (!bwi.blacklist.contains(bwi.nextRecord))
        {
            if (contains(bwi.intervals, bwi.nextRecord.rNextId, bwi.nextRecord.pNext))
            {
                unsigned rg = getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG);
                auto it = bwi.translocationBuffers[rg].pairs.find(bwi.nextRecord.qName);
                if (it != bwi.translocationBuffers[rg].pairs.end())
                    processTranslocatedRecord(bwi.translocationBuffers[rg], it, bwi.nextRecord);
            }
        }
    }
    else    // The mate cannot be in the buffer yet, so insert a new entry.
    {
        if (!bwi.blacklist.blacklisted(bwi.nextRecord))
        {
            if (contains(bwi.intervals, bwi.nextRecord.rNextId, bwi.nextRecord.pNext))
            {
                unsigned rg = getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG);
                processTranslocatedRecord(bwi.translocationBuffers[rg], bwi.nextRecord);
            }
        }
    }
}
// ---------------------------------------------------------------------------------------
// Function processRecord()
// ---------------------------------------------------------------------------------------
// Return 1 if processing of read (nextRecord) is done (=Checking requirements and maybe adding its info to bwi).
// Return 0 if a good rev. read with a good mate was found whose processing needs to wait until window is reached.
// Return -1 if the record lies outside of the interval specified in bwi.interval.
// Note: Cases where beginPos == pNext are ignored.
inline int processRecord(BamWindowIterator & bwi)
{
    if (!meetsRequirements(bwi.nextRecord, bwi.qualReq))
    {
        if (isTranslocated(bwi.nextRecord))
            checkAndRemoveTranslocatedRecord(bwi);

        return 1;
    }
    else
    {
        if (isTranslocated(bwi.nextRecord))
        {
            if (meetsTranslocRequirements(bwi.nextRecord, bwi.qualReq.minMappingQual))
                checkAndAddTranslocatedRecord(bwi);
            else
                checkAndRemoveTranslocatedRecord(bwi);

            return 1;
        }
        else if (isLeftAlignment(bwi.nextRecord))
        {
            if (processLeftAlignment(bwi))
                return 1;                           // Continue reading records.
            else
                return -1;                          // Inveral has been completely processed
        }
        else if (isRightAlignment(bwi.nextRecord))  // Not really necessary but avoids lookups when beginPos == pNext.
        {
            if (processRightAlignment(bwi))
                return 0;                           // Good pair has been found and can be written to output.
        }
        return 1;                                   // Continue reading records.
    }
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
// Find the next reverse record in the intervals and process it using processRecord. Called by goToFirstReverseRecord.
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
            else if (pr < 0 && bwi.goodLeftReads.empty())        // Record is outside of interval.
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
            int32_t chrom = (*bwi.currentReadingInterval).rID;
            ++bwi.currentReadingInterval;
            if (bwi.currentReadingInterval >= end(bwi.intervals))
            {
                return false;                               // Last interval has been processed.
            }
            if ((*bwi.currentReadingInterval).rID != chrom)
            {
                // std::cout << "rID " << chrom << " has been marked as completed after setting reading interval to "
                //           << (*bwi.currentReadingInterval).seqName << "(" << (*bwi.currentReadingInterval).rID << ")"
                //           << ":" << (*bwi.currentReadingInterval).beginPos << std::endl;
                bwi.finishedRIDs[chrom] = true;
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
    std::stable_sort(begin(bwi.intervals), end(bwi.intervals), &lowerRIDGenomicRegion);
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
    else if (pr < 0 && bwi.goodLeftReads.empty())            // Record is outside of interval.
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
// Add info about the number of overlapped windows and the tip distances of the record to the window of read group.
inline void addRecordToWindows(BamWindowIterator & bwi)
{
    //Note: bwi.nextRecords is the rightmost read of a pair at this point!
    SEQAN_ASSERT_LEQ(posToWindow(bwi.leftTipPos, bwi)->beginPos, bwi.leftTipPos);
    SEQAN_ASSERT_GT(posToWindow(bwi.leftTipPos, bwi)->beginPos + 256, bwi.leftTipPos);
    addRecord(*posToWindow(bwi.leftTipPos, bwi),
              bwi.leftTipPos,
              bwi.distance,
              bwi.clipping,
              bwi.orientation,
              getReadGroup(bwi.nextRecord.tags, bwi.readGroups, bwi.mergeRG),
              length(bwi.readGroups));
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
    int32_t maxTipDistance = max(bwi.maxTipDistances);
    int32_t tipPos = getTip(bwi.nextRecord);
    uint32_t pos = tipPos > maxTipDistance ? tipPos - maxTipDistance : 0;
    bwi.currentWindowBegin = (pos / bwi.windowShift) * bwi.windowShift;
    // Clear all windows and set their chromosome and begin position.
    for (int32_t i = bwi.currentWindowBegin;
         i < bwi.currentWindowBegin + maxTipDistance + bwi.windowSize + 1;
         i += bwi.windowShift)
        *posToWindow(i, bwi) = Window(bwi.currentWindowInterval->rID, i);
}
// ---------------------------------------------------------------------------------------
// Function isLastWindow()
// ---------------------------------------------------------------------------------------
// Return true if the last window of the (globally) last interval has been processed, false otherwise.
inline bool isLastWindow(const BamWindowIterator & bwi)
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
inline bool isLastWindowOnSeq(const BamWindowIterator & bwi)
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
        int maxTipDistance = max(bwi.maxTipDistances);
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
        //int maxTipDistances = bwi.maxTipDistances[getReadGroup(bwi.nextRecord.tags, bwi.readGroups)];
        while (bwi.currentReadingInterval == bwi.currentWindowInterval && inWindowRange(bwi, maxTipDistance))
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
// Function calculateMaxTipDistances()
// =======================================================================================
// Sets the maxTipDistances for each read group.
// Return the biggest maxTipDistance of all read groups.
inline void calculateMaxTipDistances (BamWindowIterator & bwi, const String<Histogram> & hists)
{
    unsigned rgCount = length(hists);
    resize(bwi.maxTipDistances, rgCount);
    for (unsigned rg = 0; rg < rgCount; ++rg)
    {
        bwi.maxTipDistances[rg] = hists[rg].median
                                  + 2 * hists[rg].stddev
                                  + bwi.maxDeletionSize;
    }
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
                const unsigned numContigs,
                const bool mergeRG)
{
    // Set the parameters.
    bwi.mergeRG = mergeRG;
    bwi.windowSize = windowSize;
    bwi.windowShift = windowShift;
    if (mergeRG)
    {
        bwi.readGroups[readGroups.begin()->first] = readGroups.begin()->second;
        resize(bwi.translocationBuffers, 1, Exact());
    }
    else
    {
        bwi.readGroups = readGroups;
        resize(bwi.translocationBuffers, length(hists), Exact());
    }
    // Assign readLength to BWI.
    resize(bwi.readLengths, length(hists), Exact());
    typedef Iterator<const String<Histogram>, Standard>::Type tHistIter;
    Iterator<String<int32_t> >::Type lIt = begin(bwi.readLengths, Standard());
    for (tHistIter hIt = begin(hists, Standard()); hIt != end(hists, Standard()); ++hIt, ++lIt)
        *lIt = hIt->readLength;

    bwi.maxDeletionSize = maxDeletionSize;
    calculateMaxTipDistances (bwi, hists);
    bwi.intervals = intervals;
    initContigStatusList(bwi, numContigs);
    // Set the quality requirements for reads.
    bwi.qualReq = qualReq;
    // Open the input bam file, read the header and go to the reverse read of the first good read pair in infile.
    if (!goToFirstReverseRecord(bwi))
        SEQAN_THROW(IOError("[PopDel] No BAM records meet the requirements."));
    // Initialize windows.
    int maxTipDistance = max(bwi.maxTipDistances);    //TODO: Check if this works out and replace max with parameter.
    //int maxTipDistance = bwi.maxTipDistances[getReadGroup(bwi.nextRecord.tags, bwi.readGroups)];
    resize(bwi.activeWindows, ((maxTipDistance + bwi.windowSize) / bwi.windowShift) + 1);
    resetWindows(bwi);

    // Add all records that may overlap the first window.
    while (bwi.currentReadingInterval == bwi.currentWindowInterval && inWindowRange(bwi, maxTipDistance))
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
