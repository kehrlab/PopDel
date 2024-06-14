#ifndef PROFILE_TRANSLOCATION_POPDEL_H_
#define PROFILE_TRANSLOCATION_POPDEL_H_

#include <map>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>

#include "window_popdel.h"
#include "bam_window_iterator_popdel.h"

using namespace seqan;


// =======================================================================================
// Forward declarations
// =======================================================================================
inline uint32_t getLeftClip(const BamAlignmentRecord & rec);
inline uint32_t getRightClip(const BamAlignmentRecord & rec);

// =======================================================================================
// Class TranslocationRead
// =======================================================================================
// Holds the info defining a single read of a read pair affected by a translocation
struct TranslocationRead
{
    uint32_t        refID;          // ID of the chromsome the read is mapped to
    uint32_t        pos;            // Position on the reference
    unsigned char   clip_0;           // Number of soft clipped bases at the 5'-end of the read.
    unsigned char   clip_1;           // Number of soft clipped bases at the 3'-end of the read.

    TranslocationRead(): refID(maxValue<uint32_t>()), pos(0), clip_0(0), clip_1(0) {}

    TranslocationRead(unsigned r, unsigned p, unsigned c_0, unsigned c_1): refID(r), pos(p), clip_0(c_0), clip_1(c_1) {}

    TranslocationRead& operator= (TranslocationRead other)
    {
        std::swap(refID, other.refID);
        std::swap(pos, other.pos);
        std::swap(clip_0, other.clip_0);
        std::swap(clip_1, other.clip_1);
        return *this;
    }
};
// =======================================================================================
// Class TranslocationPair
// =======================================================================================
// Holds a two reads making up a read pair affected by a translocation. Also holds the orientation.
struct TranslocationPair
{
    TranslocationRead first;
    TranslocationRead second;
    Orientation orientation;            // Orientation of first and second, in this order.

    TranslocationPair(unsigned r, unsigned p, unsigned char c_0, unsigned char c_1, Orientation o)
    {
        first.refID = r;
        first.pos = p;
        first.clip_0 = c_0;
        first.clip_1 = c_1;
        orientation = o;
    }

    TranslocationPair(const TranslocationRead & f, const TranslocationRead & s, Orientation o)
    {
        first = f;
        second = s;
        orientation = o;
    }
    inline void print()
    {
        std::cout << first.refID << ":" << first.pos << "(" << static_cast<int>(first.clip_0) << "," << static_cast<int>(first.clip_1) << ")+"
                  << second.refID << ":" << second.pos << "(" << static_cast<int>(second.clip_0) << "," << static_cast<int>(second.clip_1) << ")";
    }
};

inline bool operator<(const TranslocationPair & l, const TranslocationPair & r)
{
    if (l.first.refID > r.first.refID)
        return false;
    else if (l.first.refID < r.first.refID)
        return true;
    else if (l.first.pos  > r.first.pos)
        return false;
    else if (l.first.pos < r.first.pos)
        return true;
    else if (l.second.refID > r.second.refID)
        return false;
    else if (l.second.refID < r.second.refID)
        return true;
    else if (l.second.pos > r.second.pos)
        return false;
    else if (l.second.pos < r.second.pos)
        return true;
    else
        return l.orientation < r.orientation;
}
// =======================================================================================
// Function isOrphan()
// =======================================================================================
// Return true, if the pair.second has never been updated, i.e. pair.first is an orphan.
// Return false otherwise.
inline bool isOrphan(const TranslocationPair & pair)
{
    return pair.second.refID == maxValue<uint32_t>();
}
// Comparison function for position ranges and single positions
inline bool comesAfterRange(const Pair<unsigned> & e, const unsigned & v)
{
    return (e.i2 <= v);
}
// Comparison function for pairs of unsigned
inline bool smallerPair(const Pair<unsigned> & l, const Pair<unsigned> & r)
{
    if (l.i1 < r.i1)
        return true;
    else if (l.i1 > r.i1)
        return false;
    else
        return l.i2 < r.i2;
}
// =======================================================================================
// Class leftExtendAndMergePairs
// =======================================================================================
// Takes a sorted list of Intervals [start, end[ and merges them, if they overlap after extending them l bp to the left.
inline void leftExtendAndMergePairs(String<Pair<unsigned> > & pairs, const unsigned l)
{
    if (empty(pairs))
        return;

    unsigned i = 0;
    unsigned j = 0;
    pairs[j].i1  = pairs[j].i1 > l ? pairs[j].i1 - l : 0;
    ++j;
    while (j < length(pairs))
    {
        pairs[j].i1  = pairs[j].i1 > l ? pairs[j].i1 - l : 0;
        if (pairs[i].i2 >= pairs[j].i1)
        {
            pairs[i].i2 = pairs[j].i2;
        }
        else
        {
            ++i;
            if (i != j)
                pairs[i] = pairs[j];
        }
        ++j;
    }
    resize(pairs, i + 1, Exact());
}
// =======================================================================================
// Class TranslocationBlacklist
// =======================================================================================
// Manages the blacklist for translocations
struct TranslocationBlacklist
{
    String<String<Pair<unsigned> > > regions; // Holds a sorted String of [start,end] for each contig. One per rID

    // =======================================================================================
    // Function addBedRecord()
    // =======================================================================================
    // Take a BED record and add it to the blacklist. Assumes that 
    inline void addBedRecord(const BedRecord<Bed3> & record, const std::map<CharString, unsigned> & contigNameMap)
    {
        auto it = contigNameMap.find(record.ref);
        if (it != contigNameMap.end())     //Entries not in the contig name list of the sample can be ignored.
        {
            SEQAN_ASSERT_GEQ(length(regions), it->second);
            appendValue(regions[it->second], Pair<unsigned>(record.beginPos, record.endPos));
        }
    }
    // =======================================================================================
    // Function generateContigNameMap()
    // =======================================================================================
    // Generate an std::map mapping the contig names to their rID, i.e thei position in contigNames
    inline void generateContigNameMap(std::map<CharString, unsigned> & contigNameMap,
                                      const String<CharString> & contigNames)
    {
        for (unsigned i = 0; i < length(contigNames); ++i)
            contigNameMap[contigNames[i]] = i;
    }
    // =======================================================================================
    // Function loadFromBED()
    // =======================================================================================
    // Read the regions from the given BED file and put the entries into the blacklist.
    // If the file name is the empty string, this.regions will only be resized but left empty.
    // ReadLenght is used to extend all regions to the left by ReadLenght bp and merge intervals that overlap after
    // the extension. ReadLength should correspong to the maximum read length of all histograms, like given
    // by getMaxHistReadLen().
    inline void loadFromBED(const CharString & filename,
                            const String<CharString> & contigNames,
                            const unsigned ReadLenght)
    {
        resize(regions, length(contigNames), Exact());
        if (filename == "")
            return;

        BedFileIn bedIn;
        if (!open(bedIn,toCString(filename)))
        {
            std::ostringstream msg;
            msg << "[PopDel] Could not open BED file \'" << filename << "\'. Terminating";
            SEQAN_THROW(IOError(toCString(msg.str())));
        }
        std::map<CharString, unsigned> contigNameMap;
        generateContigNameMap(contigNameMap, contigNames);
        BedRecord<Bed3> record;
        while(!atEnd(bedIn))
        {
            try
            {
                readRecord(record, bedIn);
            }
            catch (Exception const & e)
            {
                std::ostringstream msg;
                msg << "[PopDel] Error while trying to read record from BED file \'"
                    << filename << "\':" << e.what() << ". Terminating.";
                SEQAN_THROW(IOError(toCString(msg.str())));
            }
            addBedRecord(record, contigNameMap);
        }
        unsigned count = 0;
        for (unsigned i = 0; i < length(regions); ++i)
        {
            std::sort(begin(regions[i]), end(regions[i]), smallerPair);
            leftExtendAndMergePairs(regions[i], ReadLenght);
            count += length(regions[i]);
        }
        std::ostringstream msg;
        msg << "Generated translocation blacklist with " << count
            << " regions from BED file \'" << filename << "\'.";
        printStatus(msg);
    }
    // =======================================================================================
    // Function contains()
    // =======================================================================================
    // Return true of the position given by contig and pos is contained in the blacklist, false otherwise.
    inline bool contains(const unsigned rID, const unsigned pos) const
    {
        if (length(regions) <= rID) // make sure that contig number is in range
            return false;
        auto it = std::lower_bound(begin(regions[rID]), end(regions[rID]), pos, comesAfterRange);
        if (it == end(regions[rID]))
            return false;
        else
            return pos >= it->i1 && pos < it->i2;
    }
    // =======================================================================================
    // Function contains()
    // =======================================================================================
    // Return true of the position given by contig and pos is contained in the blacklist, false otherwise.
    inline bool contains(const BamAlignmentRecord & record)
    {
        return contains(record.rID, record.beginPos);
    }
    // Returns true if one or both reads map into a blacklisted region.
    inline bool blacklisted(const BamAlignmentRecord & record)
    {
        if (contains(record.rID, record.beginPos))
            return true;
        else
            return contains(record.rNextId, record.pNext);
    }
};
inline bool comesAfterGenomicRegion(const GenomicRegion & l , const Pair<int> & r)
{
    return l.rID < r.i1 || (l.rID == r.i1 && l.beginPos < r.i2 && l.endPos <= r.i2);
}
// =======================================================================================
// Function contains()
// =======================================================================================
// Return true of the position given by contig and pos is contained in the the string of intervals. false otherwise.
inline bool contains(const String<GenomicRegion> & intervals,
                     const int32_t rID,
                     const int32_t pos)
{
    auto it = std::lower_bound(begin(intervals), end(intervals), Pair<int>(rID, pos), comesAfterGenomicRegion);
    if (it == end(intervals))
        return false;
    else
        return it->rID == rID && pos >= it->beginPos && pos < it->endPos;
}
// =======================================================================================
// Class TranslocationBuffer
// =======================================================================================
// Manages a map of Pairs of TranslocationReads
struct TranslocationBuffer
{
    std::map<CharString, TranslocationPair>      pairs;  // Maps the read ID to a read pair.

    // =======================================================================================
    // Function insert()
    // =======================================================================================
    // Insert a pair with the given values as the first read in the pair.
    inline void insert(const CharString & qname, uint32_t refID, uint32_t pos, uint32_t clip_0, uint32_t clip_1, bool reverse)
    {
        SEQAN_ASSERT_LEQ(clip_0, 255u);
        SEQAN_ASSERT_LEQ(clip_1, 255u);
        Orientation o = reverse ? Orientation::RF : Orientation::FF; // Orientation of mate is preliminary
        pairs.emplace(qname, TranslocationPair(refID, pos, clip_0, clip_1, o));
    }
    // =======================================================================================
    // Function update()
    // =======================================================================================
    // Update a pair by setting the given values for the second red in the pair.
    inline void update(std::map<CharString, TranslocationPair>::iterator it,
                       uint32_t refID,
                       uint32_t pos,
                       uint32_t clip_0,
                       uint32_t clip_1,
                       bool reverse)
    {
        SEQAN_ASSERT_LEQ(clip_0, 255u);
        SEQAN_ASSERT_LEQ(clip_1, 255u);
        it->second.second = TranslocationRead(refID, pos, clip_0, clip_1);
        if (it->second.orientation == Orientation::FF)  // Must match orientation in insert()!
            it->second.orientation = reverse ? Orientation::FR : Orientation::FF;
        else
            it->second.orientation = reverse ? Orientation::RR : Orientation::RF;
    }
    // =======================================================================================
    // Function removePair()
    // =======================================================================================
    // Remove the pair the iterator points at.
    // Update the iterator to point at the next element after the deleted one.
    template<typename TIter>
    inline void removePair(TIter& it)
    {
        it = pairs.erase(it);
    }
    // =======================================================================================
    // Function removeIfFound()
    // =======================================================================================
    // Check if a pair with given rID already exists in the map.
    // If no, do nothing.
    // If yes, remove it.
    inline void removeIfFound(const CharString & qname)
    {
        auto it = pairs.find(qname);
        if (it != pairs.end())
            removePair(it);
    }
    // =======================================================================================
    // Function purgeOrphans()
    // =======================================================================================
    // Remove all orphans, i.e. "pairs" without an updated second entry.
    inline void purgeOrphans()
    {
        for (std::map<CharString, TranslocationPair>::const_iterator it = pairs.begin(); it != pairs.end();)
            if (isOrphan(it->second))
                removePair(it);
            else
                ++it;
    }
    inline void report()
    {
        std::cout << "Translocation Buffer report" << std::endl;
        std::cout << "#Pairs:\t" << pairs.size() << std::endl;
        for (auto it = pairs.begin(); it != pairs.end(); ++it)
        {
            std::cout << it->first << "\t";
            it->second.print();
            std::cout << std::endl;
        }
    }
};
// =======================================================================================
// Function moveTranslocationBufferToString()
// =======================================================================================
// Moves the entries of the TranslocationBuffer source s to the target string t.
// Also creates the reverse entry for each translocation pair.
inline void moveTranslocationBufferToString(String<TranslocationPair> & t, TranslocationBuffer & s)
{
    reserve(t, 2 * s.pairs.size(), Exact());
    auto it = s.pairs.begin();
    while (it != s.pairs.end())
    {
        appendValue(t, it->second);
        appendValue(t, TranslocationPair(it->second.second,
                                         it->second.first,
                                         getSwitchedOrientation(it->second.orientation)));
        it = s.pairs.erase(it);
    }
}
// =======================================================================================
// Function isTranslocated()
// =======================================================================================
// Return true if the record's other read has been mapped to another chromsome.
inline bool isTranslocated(const BamAlignmentRecord & record)
{
    return record.rID != record.rNextId;
}
// =======================================================================================
// Function meetsTranslocRequirements()
// =======================================================================================
// Return true if the record meets the additional requirements for translocation reads.
// Return false otherwise.
inline bool meetsTranslocRequirements(const BamAlignmentRecord & record, const unsigned minMappingQual = 30)
{
    return record.mapQ >= minMappingQual;
}
// =======================================================================================
// Function createTranslocationBufferEntry()
// =======================================================================================
// Determine the parameters of a translocationBuffer entry from the record.
inline void createTranslocationBufferEntry(const BamAlignmentRecord & record,
                                           unsigned & pos,
                                           unsigned & clip_0,
                                           unsigned & clip_1,
                                           bool & reverse)
{
    reverse = hasFlagRC(record);
    if (reverse)
    {
        pos = record.beginPos;
        clip_0 = getRightClip(record);
        clip_1 = getLeftClip(record);
    }
    else
    {
        pos = record.beginPos + getAlignmentLengthInRef(record) - 1;
        clip_0 = getLeftClip(record);
        clip_1 = getRightClip(record);
    }
}
// =======================================================================================
// Function processTranslocatedRecord()
// =======================================================================================
// Processes a record whose other read has been mapped to another chromosome and adds it to the buffer.
inline void processTranslocatedRecord(TranslocationBuffer & buffer, const BamAlignmentRecord & record)
{
    unsigned pos;
    unsigned clip_0;
    unsigned clip_1;
    bool reverse;
    createTranslocationBufferEntry(record, pos, clip_0, clip_1, reverse);
    buffer.insert(record.qName, record.rID, pos, clip_0, clip_1, reverse);
}
// Overload for updating the record instead of inserting a new one
inline void processTranslocatedRecord(TranslocationBuffer & buffer,
                                      std::map<CharString, TranslocationPair>::iterator it,
                                      const BamAlignmentRecord & record)
{
    unsigned pos;
    unsigned clip_0;
    unsigned clip_1;
    bool reverse;
    createTranslocationBufferEntry(record, pos, clip_0, clip_1, reverse);
    buffer.update(it, record.rID, pos, clip_0, clip_1, reverse);
}
// =======================================================================================
// Function getTranslocationWindowBorder()
// =======================================================================================
// Return an iterator pointing to the first element between start and end whose pos no longer within winStartPos + 256 bp.
// Count the number of records in counter.
typedef typename Iterator<String<TranslocationPair> >::Type TTranSIt;
inline TTranSIt getTranslocationWindowBorder(uint32_t & counter,
                                             const TTranSIt start,
                                             const TTranSIt end,
                                             const uint32_t chrom,
                                             const uint32_t winStartPos)
{
    TTranSIt it = start;
    const uint32_t border = winStartPos + 256;
    counter = 0;
    while (it != end && it->first.refID == chrom && it->first.pos < border)
    {
        ++it;
        ++counter;
    }
    return it;
}
// =======================================================================================
// Function packOrientations()
// =======================================================================================
// Overload of packOrientations() in window_popdel.h. Uses for tranlocation strings.
// Put the orientations in groups of 4, s.t. they can fill all 8 bits of a uint8_t. Put these uints8_t into an a String.
// Return the length of the resulting string.
inline unsigned packOrientations(String<uint8_t> & packed,
                                 const TTranSIt start,
                                 const TTranSIt end)
{
    unsigned numRecords = std::distance(start, end);
    unsigned numChars = std::ceil(static_cast<double>(numRecords) / 4);
    resize(packed, numChars, 0, Exact());
    unsigned i = 0;
    for (auto it = start; it != end; ++it)
    {
        unsigned j = i / 4;    // Fill all 8 bits of the current uint8_t before moving to the next one.
        packed[j] <<= 2;       // Create space for the 2 new orientation bits.
        packed[j] |= static_cast<uint8_t>(it->orientation);
        ++i;
    }
    packed[numChars-1] <<= numChars * 8 - numRecords * 2;    // Catch up on the remaining shifts, if there are any.
    return numChars;
}
// =======================================================================================
// Function prepareTranslocationBlock()
// =======================================================================================
inline void prepareTranslocationBlock(std::ofstream & out,
                                      String<String<TranslocationPair> > & bufferStrings,
                                      String<TranslocationBuffer> & buffers)
{
    // std::cout << "Translocation buffer is holding records of " << length(buffers) << " RG(s):" << std::endl;
    for (unsigned rg = 0; rg < length(buffers); ++rg)
    {
        // std::cout << "Size of translocation map of RG " << rg << " before removing orphans: " << buffers[rg].pairs.size() << " entries." << std::endl;
        buffers[rg].purgeOrphans();
        // std::cout << "Size of translocation map of RG " << rg << " after removing orphans: " << buffers[rg].pairs.size() << " entries:" << std::endl;
    }
    resize(bufferStrings, length(buffers));
    for (unsigned i = 0; i < length(buffers); ++i)
    {
        moveTranslocationBufferToString(bufferStrings[i], buffers[i]);
        std::sort(begin(bufferStrings[i]), end(bufferStrings[i]));
    }
    resize(buffers, 0);

    uint32_t guard = maxValue<uint32_t>();           // Must be of the same type as 'chrom' in writeWindow()!
    zlib_stream::zip_ostream stream(out);
    stream.write(reinterpret_cast<char *>(&guard), sizeof(uint32_t));
    stream.zflush();
}
// =======================================================================================
// Function getMinStartingPos()
// =======================================================================================
// Get the first starting position among all read groups and prepare the itString.
inline void getMinStartingPos(String<TTranSIt> & itString,
                              String<String<TranslocationPair> > & bufferStrings,
                              uint32_t & chrom,
                              uint32_t & pos)
{
    for (unsigned rg = 0; rg < length(itString); ++rg)
    {
        itString[rg] = begin(bufferStrings[rg]);
        if (empty(bufferStrings[rg]))
            continue;
        if (chrom > itString[rg]->first.refID)
        {
            chrom = itString[rg]->first.refID;
            pos = itString[rg]->first.pos;
        }
        else if (chrom == itString[rg]->first.refID && pos > itString[rg]->first.pos)
        {
            pos = itString[rg]->first.pos;
        }
    }
}
// =======================================================================================
// Function writeWindowRecords()
// =======================================================================================
// Write the translocation records for one RG between first and last for the given winStartPos to stream.
// Also move first
inline void writeWindowRecords(zlib_stream::zip_ostream & stream,
                               TTranSIt & first,
                               const TTranSIt last,
                               const unsigned winStartPos)
{
    for (;first != last; ++first)
    {
        char offset = first->first.pos - winStartPos;
        unsigned char clip_0 = first->first.clip_0;
        unsigned char clip_1 = first->first.clip_1;
        uint32_t refID2 = first->second.refID;
        uint32_t pos2 = first->second.pos;
        unsigned char clip_2 = first->second.clip_0;
        unsigned char clip_3 = first->second.clip_1;
        stream.write(reinterpret_cast<char *>(&offset), sizeof(char));
        stream.write(reinterpret_cast<char *>(&clip_0), sizeof(unsigned char));
        stream.write(reinterpret_cast<char *>(&clip_1), sizeof(unsigned char));
        stream.write(reinterpret_cast<char *>(&refID2), sizeof(uint32_t));
        stream.write(reinterpret_cast<char *>(&pos2), sizeof(uint32_t));
        stream.write(reinterpret_cast<char *>(&clip_2), sizeof(unsigned char));
        stream.write(reinterpret_cast<char *>(&clip_3), sizeof(unsigned char));
    }
}
// =======================================================================================
// Function maybeUpdateTranslocationIndexFields()
// =======================================================================================
// Update the TranslocationIndexFields if the chromsome has changed.
// Update refID1 if necessary.
inline void maybeUpdateTranslocationIndexFields(std::ofstream & out,
                                                String<uint64_t> & translocationIndexFields,
                                                uint32_t & refID1,
                                                const uint32_t chrom)
{
    if (refID1 != chrom)
    {
        refID1 = chrom;
        //std::cout << "Writing file position " << out.tellp() << " to index for refID " << refID1 << std::endl;
        translocationIndexFields[refID1] = out.tellp();
        //std::cout << "." << std::endl;
    }
}
// =======================================================================================
// Function updateChromAndPos()
// =======================================================================================
// Return true if there are more windows to process, false otherwise.
// Also update chrom and firstPos accordingly.
inline bool updateChromAndPos(uint32_t & chrom,
                              uint32_t & firstPos,
                              const TTranSIt it,
                              const TTranSIt ending)
{
    if (it != ending)
    {
        if (it->first.refID < chrom)
        {
            chrom = it->first.refID;
            firstPos = it->first.pos;
            return true;
        }
        else if (it->first.refID == chrom && it->first.pos < firstPos)
        {
             firstPos =  it->first.pos;
             return true;
        }
    }
    return false;
}
// =======================================================================================
// Function writeTranslocationBlock()
// =======================================================================================
inline void writeTranslocationBlock(std::ofstream & out,
                                    String<uint64_t> & translocationIndexFields,
                                    String<TranslocationBuffer> & buffers)
{
    // for (unsigned i = 0; i < length(buffers); ++i)
    // {
    //     buffers[i].report();
    // }
    String<String<TranslocationPair> > bufferStrings;
    prepareTranslocationBlock(out, bufferStrings, buffers); // TODO Continue here and check if invalid entries are in the Block
    String<TTranSIt> itString;
    resize(itString, length(bufferStrings), Exact());
    uint32_t chrom = maxValue<uint32_t>();
    uint32_t firstPos = maxValue<uint32_t>();
    // Get the min starting position of the first record of all read groups.
    getMinStartingPos(itString, bufferStrings, chrom, firstPos);
    bool moreWindowsToProcess = chrom != maxValue<uint32_t>();
    uint32_t refID1 = maxValue<uint32_t>();
    while (moreWindowsToProcess)    // While there are still windows to process ...
    {
        maybeUpdateTranslocationIndexFields(out, translocationIndexFields, refID1, chrom);
        zlib_stream::zip_ostream stream(out);
        // Write Chromsome and Position of Window
        uint32_t winStartPos = (firstPos / 256) * 256;
        stream.write(reinterpret_cast<char *>(&refID1), sizeof(uint32_t));
        stream.write(reinterpret_cast<char *>(&winStartPos), sizeof(uint32_t));
        firstPos = maxValue<unsigned>();    // Set to maxValue in preparation of getting the next window pos.
        chrom = maxValue<unsigned>();
        moreWindowsToProcess = false;   // Will become true if getTranslocationWindowBorder() does not return end() for all RG.
        for (unsigned rg = 0; rg < length(bufferStrings); ++rg)
        {
            uint32_t numReadPairs;
            // Get the border of the current window by moving itString[rg] to the first entry behind it.
            TTranSIt it = itString[rg];
            TTranSIt ending = end(bufferStrings[rg]);
            itString[rg] = getTranslocationWindowBorder(numReadPairs, it, ending, refID1, winStartPos);
            // Write the number of read pairs in this read group
            stream.write(reinterpret_cast<char *>(&numReadPairs), sizeof(uint32_t));
            if (numReadPairs != 0)
            {
                // Pack the orientations and write them
                String<uint8_t> orientations;
                unsigned numChars = packOrientations(orientations, it, itString[rg]);
                stream.write(reinterpret_cast<char *>(&orientations[0]), numChars);
                writeWindowRecords(stream, it, itString[rg], winStartPos);
            }
            // Check if there are still more windows to process
            moreWindowsToProcess |= updateChromAndPos(chrom, firstPos, it, ending);
        }
        stream.zflush();
    }
}
// =======================================================================================
// Function writeTranslocationIndexIntoHeader()
// =======================================================================================
// Write the fields of the translocation index to the stream. Replace empty index fields with maxOffset
// Note: Assumes that the stream is already pointing to the start position of the translocation index.
template<typename TStream>
inline void writeTranslocationIndexIntoHeader(TStream & stream,
                                              String<uint64_t> translocationIndex,
                                              uint64_t maxOffset)
{
    // Fill in file offsets for empty regions.
    uint64_t prev = maxOffset;
    for (int i = length(translocationIndex) - 1; i >= 0; --i)
    {
        if (translocationIndex[i] == 0)
            translocationIndex[i] = prev;
        else
            prev = translocationIndex[i];
    }
    for (unsigned i = 0; i < length(translocationIndex); ++i)
        stream.write(reinterpret_cast<char *>(&translocationIndex[i]), sizeof(uint64_t));
}
#endif /* PROFILE_TRANSLOCATION_POPDEL_H_ */