#ifndef WINDOW_POPDEL_H_
#define WINDOW_POPDEL_H_

#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include "../histogram_popdel.h"

using namespace seqan;

// =======================================================================================
// Enum Struct Orientation
// =======================================================================================
// Enum for orientation of the read pairs.
enum struct Orientation: uint8_t
{
    FR = 0,     // > <
    FF = 1,     // > >
    RR = 2,     // < <
    RF = 3      // < >
};
// =======================================================================================
// Function firstReadForward()
// =======================================================================================
// Return true if the first read of the read pair is pointing forward.
// Return false otherwise
inline bool firstReadForward(const Orientation & o)
{
    return o == Orientation::FR || o == Orientation::FF;
}
// =======================================================================================
// Function secondReadForward()
// =======================================================================================
// Return true if the second read of the read pair is pointing forward.
// Return false otherwise
inline bool secondReadForward(const Orientation & o)
{
    return o == Orientation::RF || o == Orientation::FF;
}
// =======================================================================================
// Function getPairOrientation()
// =======================================================================================
// Return the orientation of the read pair.
// IMPORTANT: Assumes that the given record is the rightmost read of the pair!
inline Orientation getPairOrientation(const uint32_t flag)
{
    if (flag & 16)           // read reverse strand
    {
        if (flag & 32)       // mate reverse strand
        {
            return Orientation::RR;
        }
        else
        {
            return  Orientation::FR;
        }
    }
    else if (flag & 32)      // mate reverse strand
    {
        return Orientation::RF;
    }
    else
    {
        return Orientation::FF;
    }
}
// Overload for directly applying getPairOrientation() on a BamAlignmentRecord.
inline Orientation getPairOrientation(const BamAlignmentRecord & record)
{
    return getPairOrientation(record.flag);
}
std::ostream & operator << (std::ostream & out, const Orientation & o)
{
    if (o == Orientation::FR )
        out << "FR";
    else if (o == Orientation::FF )
        out << "FF";
    else if (o == Orientation::RR )
        out << "RR";
    else
        out << "RF";
    return out;
}
// =======================================================================================
// Function getSwitchedOrientation()
// =======================================================================================
// Return the orientation one would get by switching first and second read, but keeping their individual orientations.
// FR -> RF, RF -> FR, but FF -> FF and RR -> RR.
inline Orientation getSwitchedOrientation(const Orientation o)
{
    if (o == Orientation::FR)
        return Orientation::RF;
    else if (o == Orientation::RF)
        return Orientation::FR;
    else
        return o;
}
// =======================================================================================
// Struct ReadPair
// =======================================================================================
// Holds the properties of one record (=read pair) in the window.
struct ReadPair
{
    uint32_t pos;                   // Position of the 3' end of the leftmost read in the pair
    int32_t distance;                   // FR-reads: 5'-end distance. Other reads; 3'-end distance
    std::vector<uint8_t> clipping;   // number of clipped bases (5', left; 3', left; 3', right; 5', right)
    uint16_t totalClipping;
    Orientation orientation;

    ReadPair():
        pos(0), distance(0), clipping(std::vector<uint8_t>(4, 0)), totalClipping(0), orientation(Orientation::FR )
    {}

    ReadPair(unsigned p, unsigned i, uint8_t c, Orientation o = Orientation::FR ):
        pos(p), distance(i), clipping({0, c/2, c/2, 0}), totalClipping(c), orientation(o)
    {}

    ReadPair(unsigned p, unsigned i, uint8_t c_0, uint8_t c_1, uint8_t c_2, uint8_t c_3, Orientation o = Orientation::FR ):
        pos(p), distance(i), clipping({c_0, c_1, c_2, c_3}), totalClipping(c_0 + c_1 + c_2 + c_3), orientation(o)
    {}

    ReadPair(unsigned p, unsigned i, std::vector<uint8_t> c, Orientation o = Orientation::FR ):
        pos(p), distance(i), clipping(c), totalClipping(0), orientation(o)
    {
        for (unsigned j = 0; j < clipping.size(); ++j)
            totalClipping += clipping[j];
    }

    bool operator <(const ReadPair & r)
    {
        if (pos < r.pos)
            return true;
        if (pos == r.pos && distance < r.distance )
            return true;
        else
            return false;
    }
};
inline bool operator==(const ReadPair& l, const ReadPair& r)
{
    return l.pos == r.pos && l.distance == r.distance && l.clipping[0] == r.clipping[0] && l.clipping[1] == r.clipping[1] && l.clipping[2] == r.clipping[2] && l.clipping[3] == r.clipping[3] && l.orientation == r.orientation;
}
inline bool operator!=(const ReadPair& l, const ReadPair& r)
{
    return !(l == r);
}
// =======================================================================================
// Struct Window
// =======================================================================================
// Defines a window on the Genome. Can contain samples from multiple read groups.
struct Window
{
    int32_t chrom;
    int32_t beginPos;
    String<String<ReadPair> > records;

    Window():
        chrom(GenomicRegion::INVALID_ID), beginPos(GenomicRegion::INVALID_POS)
    {}

    Window(__int32 c, __int32 pos):
        chrom(c), beginPos(pos)
    {}
};

// =======================================================================================
// Function isEmpty()
// =======================================================================================
// Return true if the window is empty (i.e window.records has length 0)
template<typename TWindow>
inline bool isEmpty(const TWindow & window)
{
    return length(window.records) == 0;
}

// =======================================================================================
// Function reset()
// =======================================================================================
// Reset chrom, beginPos and records window but keep the same number of read groups.
template<typename TWindow>
inline void reset(TWindow & window)
{
    window.chrom = 0;
    window.beginPos = 0;
    for (unsigned rg = 0; rg < length(window.records); ++rg)
        clear(window.records[rg]);
}

// =======================================================================================
// Function addRecord()
// =======================================================================================
// Aggregate tipDistance and numWindows as pair and add it to window for read group.
inline void addRecord(Window & window,
                      int32_t beginPos,     // clipped 3'-end of the leftmost read in the pair
                      int32_t tipDistance,
                      std::vector<uint8_t> clipping,
                      Orientation orientation,
                      unsigned readGroup,
                      unsigned numReadGroups)
{
    SEQAN_ASSERT_NEQ(numReadGroups, 0u);
    SEQAN_ASSERT_LT(readGroup, numReadGroups);
    if (isEmpty(window))                                // TODO: Maybe move this to some initialization step.
        resize(window.records, numReadGroups);
    appendValue(window.records[readGroup], ReadPair(beginPos, tipDistance, clipping, orientation));
    SEQAN_ASSERT_LEQ(beginPos, window.beginPos + 255);
    SEQAN_ASSERT_GEQ(beginPos, window.beginPos);
    SEQAN_ASSERT_EQ(isEmpty(window), false);
}
// =======================================================================================
// Function countClippedPairs()
// =======================================================================================
// Return the number of read pairs with clipping in window for read group rg.
inline unsigned countClippedPairs(const Window & window, unsigned rg)
{
    unsigned c = 0;
    for (unsigned  i = 0; i < length(window.records[rg]); ++i)
    {
        if (window.records[rg][i].totalClipping > 0u)
            ++c;
    }
    return c;
}
// =======================================================================================
// Function countNonFRPairs()
// =======================================================================================
// Return the number of read pair in records that have a non-FR orientation.
inline unsigned countNonFRPairs(const Window & window, unsigned rg)
{
    unsigned c = 0;
    for (unsigned  i = 0; i < length(window.records[rg]); ++i)
        if (window.records[rg][i].orientation != Orientation::FR)
            ++c;
    return c;
}
// =======================================================================================
// Function packOrientations()
// =======================================================================================
// Put the orientations in groups of 4, s.t. they can fill all 8 bits of a uint8_t. Put these uints8_t into an a String.
// Ignore read pairs with forward-reverse orientation.
// Return the length of the resulting string.
inline unsigned packOrientations(String<uint8_t> & packed,
                                 String<uint32_t> & ids,
                                 const String<ReadPair> & records,
                                 const unsigned & numNonFRRecords)
{
    SEQAN_ASSERT_GT(numNonFRRecords, 0u);
    SEQAN_ASSERT(empty(ids));
    unsigned numChars = std::ceil(static_cast<double>(numNonFRRecords) / 4);
    resize(packed, numChars, 0, Generous()); // TODO: Why Generous()? Try Exact()
    reserve(ids, numNonFRRecords, Exact());
    unsigned skipped = 0;
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (records[i].orientation != Orientation::FR)
        {
            SEQAN_ASSERT_LEQ(skipped, i);
            SEQAN_ASSERT_LEQ(static_cast<unsigned>(i - skipped), length(ids));
            appendValue(ids, i);
            unsigned j = (i - skipped) / 4;     // Fill all 8 bits of the current uint8_t before moving to the next one.
            packed[j] <<= 2;                                             // Create space for the 2 new orientation bits.
            packed[j] |= static_cast<uint8_t>(records[i].orientation);
        }
        else
        {
            ++skipped;
        }
    }
    packed[numChars-1] <<= numChars * 8 - numNonFRRecords * 2;    // Catch up on the remaining shifts, if there are any.
    return numChars;
}
// =======================================================================================
// Function unpackOrientations()
// =======================================================================================
// Put the orientations in groups of 4, s.t. they can fill all 8 bits of a uint8_t. Put these uints8_t into a String.
inline void unpackOrientations(String<ReadPair> & records, const String<uint8_t> & packed, const String<unsigned> & ids)
{
    uint8_t mask = 192;     // Translates to 11000000 in binary format.
    for (unsigned i = 0; i < length(ids); ++i)
        records[ids[i]].orientation = static_cast<Orientation>((packed[i / 4] & (mask >> i % 4 * 2)) >> (6-i%4*2));
}
// =======================================================================================
// Function unpackOrientations()
// =======================================================================================
// Return the stating position of the window that contains the given position.
inline unsigned posToWin(const unsigned & pos, const unsigned windowSize = 30u)
{
    return (pos / windowSize) * windowSize;
}
// =======================================================================================
// Function writeWindow()
// =======================================================================================
// Used for writing the compressed windows. Used by PopDel profile.
inline void writeWindow(zlib_stream::zip_ostream & stream,
                        const Window & window)
{
    typedef Iterator<const String<ReadPair>, Standard>::Type TValueIter;
    SEQAN_ASSERT(!isEmpty(window));

    // Write chromosome and position.
    uint32_t chrom = window.chrom;
    stream.write(reinterpret_cast<char *>(&chrom), sizeof(uint32_t));
    uint32_t beginPos = window.beginPos;
    stream.write(reinterpret_cast<char *>(&beginPos), sizeof(uint32_t));

    for (unsigned rg = 0; rg < length(window.records); ++rg)
    {
        // Write number of read pairs for this read group.
        uint32_t numRecords = length(window.records[rg]);
        stream.write(reinterpret_cast<char *>(&numRecords), sizeof(uint32_t));

        // Write number of clipped read pairs for this read group.
        uint32_t numClippedRecords = countClippedPairs(window, rg);
        stream.write(reinterpret_cast<char *>(&numClippedRecords), sizeof(uint32_t));

        // Write number of non-FR orientated read pairs for this read group.
        uint32_t numNonFR = countNonFRPairs(window, rg);
        stream.write(reinterpret_cast<char *>(&numNonFR), sizeof(uint32_t));

        // Write the ammount of clipped bases for all read pairs that have clipping.
        // Identify them by ther position in the window
        if (numClippedRecords > 0)
        {
            uint32_t id = 0;
            for (TValueIter it = begin(window.records[rg], Standard());
                 it != end(window.records[rg], Standard());
                 ++it, ++id)
                 {
                     if (it->totalClipping > 0u)
                     {
                        stream.write(reinterpret_cast<char *>(&id), sizeof(uint32_t));
                        for (unsigned cl = 0; cl < it->clipping.size(); ++cl)
                        {
                            SEQAN_ASSERT_LEQ(it->clipping[cl], 255u);
                            unsigned char clipping = it->clipping[cl];
                            stream.write(reinterpret_cast<char *>(&clipping), sizeof(unsigned char));
                        }
                     }
                 }
        }
        // Write the packed string of read pair orientations.
        if (numNonFR > 0)
        {
            String<uint8_t> orientations;
            String<uint32_t> oIDs;
            unsigned numChars = packOrientations(orientations, oIDs, window.records[rg], numNonFR);
            stream.write(reinterpret_cast<char *>(&oIDs[0]), sizeof(oIDs[0]) * numNonFR);
            stream.write(reinterpret_cast<char *>(&orientations[0]), numChars);
        }

        // Now write the records      
        for (TValueIter it = begin(window.records[rg], Standard()); it != end(window.records[rg], Standard()); ++it)
        {
            if (static_cast<unsigned>(it->pos) > beginPos + 255u || static_cast<unsigned>(it->pos) < beginPos)
            {
                SEQAN_THROW(ParseError("Invalid offset in Window!"));
            }

            char posOffset = it->pos - window.beginPos;
            int32_t distance = it->distance;
            stream.write(reinterpret_cast<char *>(&posOffset), sizeof(char));
            stream.write(reinterpret_cast<char *>(&distance), sizeof(int32_t));
        }
    }
}
// =======================================================================================
// Function readCoordinates()
// =======================================================================================
// Read the chromosome and position of the window from the input stream. Used by readWindow().
// Return true on success, false if EOF has been reached or the translocation guard is encountered
inline bool readCoordinates(int32_t & chrom, int32_t & beginPos, zlib_stream::zip_istream & stream)
{
    // Read chromosome and position.
    uint32_t refID = chrom;
    stream.read(reinterpret_cast<char *>(&refID), sizeof(uint32_t));
    if (refID == maxValue<uint32_t>() || stream.eof())
        return false;
    SEQAN_ASSERT_LEQ(refID, static_cast<uint32_t>(maxValue<int32_t>()));
    chrom = refID;
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig name."));
    stream.read(reinterpret_cast<char *>(&beginPos), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig position."));
    return true;
}
// =======================================================================================
// Function readNumRecords()
// =======================================================================================
// Return the number of read pairs in the current window of the input stream. Used by readWindow().
inline unsigned readNumRecords(zlib_stream::zip_istream & stream)
{
    // Read the number of read pairs for this read group.
    uint32_t numRecords = 0;
    stream.read(reinterpret_cast<char *>(&numRecords), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read number of read pairs."));
    return numRecords;
}
// =======================================================================================
// Function readNumClippedRecords()
// =======================================================================================
// Return the number of clipped read pairs in the current window of the input stream. Used by readWindow().
inline unsigned readNumClippedRecords(zlib_stream::zip_istream & stream)
{
    // Read the number of read pairs for this read group.
    uint32_t numClippedRecords = 0;
    stream.read(reinterpret_cast<char *>(&numClippedRecords), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read number of clipped read pairs."));
    return numClippedRecords;
}
// =======================================================================================
// Function readNumClippedRecords()
// =======================================================================================
// Return the number of clipped read pairs in the current window of the input stream. Used by readWindow().
inline unsigned readNumNonFRRecords(zlib_stream::zip_istream & stream)
{
    // Read the number of read pairs for this read group.
    uint32_t numNonFRRecords = 0;
    stream.read(reinterpret_cast<char *>(&numNonFRRecords), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read number non-FR oriented records."));
    return numNonFRRecords;
}
// =======================================================================================
// Function readOrientations()
// =======================================================================================
// Read the packed string of read pair orientations from the stream and unpack them into the records.
// Used by readWindow().
inline void readOrientations(String<ReadPair> & records, zlib_stream::zip_istream & stream, const unsigned & numRecords)
{
    if (numRecords == 0)
        return;

    String<uint8_t> orientations;
    String<uint32_t> ids;
    unsigned orientationChars = std::ceil(static_cast<double>(numRecords) / 4);
    resize(orientations, orientationChars, Generous());
    resize(ids, numRecords, Exact());
    stream.read(reinterpret_cast<char *>(&ids[0]), sizeof(uint32_t) * numRecords);
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read orientation IDs of read pairs."));
    stream.read(reinterpret_cast<char *>(&orientations[0]), orientationChars);
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read orientations of read pairs."));
    unpackOrientations(records, orientations, ids);   // Unpack the orientations and store them in the window.
}
// =======================================================================================
// Function readClipping()
// =======================================================================================
// Read all clipping info of the window's records and create a map of ID:clipping.
inline void readClipping(std::map<unsigned, std::vector<uint8_t>> & clipMap,
                         zlib_stream::zip_istream & stream,
                         const unsigned & numClippedRecords)
{
    for (unsigned i = 0; i < numClippedRecords; ++i)
    {
        unsigned id = 0;
        std::vector<uint8_t> clipping(4, 0);

        stream.read(reinterpret_cast<char *>(&id), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read id of clipped read pair."));
        for (unsigned i = 0; i < clipping.size(); ++i)
        {
            uint8_t tempClip = 0;
            stream.read(reinterpret_cast<char *>(&tempClip), sizeof(unsigned char));
            if (!stream.good())
                SEQAN_THROW(ParseError("[PopDel] Unable to read number of clipped bases."));
            clipping[i] = tempClip;
        }
        
        clipMap[id] = clipping;
    }
}
// =======================================================================================
// Function readRecords()
// =======================================================================================
// Read all records of one read group in the current window of the input stream. Used by readWindow().
inline void readRecords(String<ReadPair> & records,
                        zlib_stream::zip_istream & stream,
                        const unsigned & windowBeginPos,
                        const unsigned & numRecords)
{
    for (uint32_t i = 0; i < numRecords; ++i)
    {
        stream.read(reinterpret_cast<char *>(&records[i].pos), sizeof(char));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read position offset."));
        records[i].pos += windowBeginPos;
        stream.read(reinterpret_cast<char *>(&records[i].distance), sizeof(int32_t));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read tip distance deviation."));
    }
}
// =======================================================================================
// Function addClippingToRecords()
// =======================================================================================
// Add the clipping info from the map to the records.
inline void addClippingToRecords(String<ReadPair> & records,
                                 const std::map<unsigned, std::vector<uint8_t>> & clipMap,
                                 const unsigned & numClippedRecords)
{
    if (numClippedRecords == 0)
        return;
    unsigned  processed = 0;
    const std::map<unsigned, std::vector<uint8_t>>::const_iterator mapEnd = clipMap.end();
    for (unsigned  i = 0; i < length(records) && processed < numClippedRecords; ++i)
    {
        std::map<unsigned, std::vector<uint8_t>>::const_iterator it = clipMap.find(i);
        if (it != mapEnd)
        {
            records[i].clipping = it->second;
            ++processed;
        }
    }
}
// =======================================================================================
// Function readWindow()
// =======================================================================================
// Read the entires of the next window in the opened and unzipped stream.
// Return true on success and false if EOF has been reached.
inline bool readWindow(Window & window,
                       zlib_stream::zip_istream & stream,
                       const unsigned & numReadGroups)
{
    reset(window);
    if (!readCoordinates(window.chrom, window.beginPos, stream))
        return false;

    resize(window.records, numReadGroups);
    for (unsigned rg = 0; rg < numReadGroups; ++rg)
    {
        std::map<unsigned, std::vector<uint8_t>> clipMap;
        unsigned numRecords = readNumRecords(stream);
        unsigned numClippedRecords = readNumClippedRecords(stream);
        unsigned numNonFRRecords = readNumNonFRRecords(stream);
        resize(window.records[rg], numRecords);
        readClipping(clipMap, stream, numClippedRecords);
        readOrientations(window.records[rg], stream, numNonFRRecords);
        readRecords(window.records[rg], stream, window.beginPos, numRecords);
        addClippingToRecords(window.records[rg], clipMap, numClippedRecords);
    }
    return true;
}
// =======================================================================================
// Function convertWindow()
// =======================================================================================
// Convert a window of size "oldWinSize" to a string of windows with "newWinSize".  //TODO: Adapt to sorted orig. win.
// Return true on success, false if there are no entries in the window.
inline bool convertWindow(String<Window> & c, const Window & o, unsigned oldWinSize = 256, unsigned newWinSize = 30)
{
    SEQAN_ASSERT_LT(newWinSize, oldWinSize);
    const unsigned rgNum = length(o.records);
    // Maxtrix used for counting how many elements will be in each window. First index: RG, second index: Window
    String<String<unsigned> > winCountMap;
    resize(winCountMap, rgNum, Exact());
    for (unsigned rg = 0; rg < rgNum; ++rg)
        resize(winCountMap[rg], std::ceil(static_cast<double>(oldWinSize + newWinSize) / newWinSize ), 0 , Exact());

    // Get the lowest starting position of all read groups
    unsigned beginPos = maxValue<unsigned>();
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        if (length(o.records[rg]) != 0)
        {
            if (o.records[rg][0].pos < beginPos)
                beginPos = o.records[rg][0].pos ;
        }
    }
    beginPos = (beginPos / newWinSize) * newWinSize;

    String<bool> windowsWithContent;
    resize(windowsWithContent, std::ceil(static_cast<double>(oldWinSize + newWinSize) / newWinSize),false, Exact());
    // Now count how many entries there will be in each window for each read group.
    // Also keep track of the number of different windows.
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        for (unsigned i = 0; i < length(o.records[rg]); ++i)
        {
            const unsigned currentWin = o.records[rg][i].pos;
            const unsigned winIndex = (currentWin - beginPos) / newWinSize;
            SEQAN_ASSERT_GEQ(winIndex, 0u);
            SEQAN_ASSERT_LEQ(winIndex, length(windowsWithContent)); // This happens iff the profiles are not sorted!
            if (winCountMap[rg][winIndex] == 0)
            {
                windowsWithContent[winIndex] = true;
            }
            ++winCountMap[rg][winIndex];
        }
    }
    // Count the number of windows with content.
    unsigned totalWinNum = 0;
    for (unsigned i = 0; i < length(windowsWithContent); ++i)
    {
        if (windowsWithContent[i])
        {
            ++totalWinNum;
        }
    }
    if (totalWinNum == 0)
    {
        return false;               // No content. Profile corrupted?
    }
    // Create a string indicating how many windows were skipped before to a given window.
    // E.g. winCountMap hold 9 windows, but 4,6,7 (starting from 0) are empty.
    // Thus cumWindowsWithContent will hold [0, 0, 0, 0, 0, 1, 1 , 2, 3]
    // We need this for adjusting the indew when saving the entries in c.
    String<unsigned> cumSkippedWindows;
    resize(cumSkippedWindows, length(windowsWithContent), 0, Exact());
    for (unsigned i = 1; i < length(cumSkippedWindows); ++i)
    {
        if (windowsWithContent[i-1])
            cumSkippedWindows[i] = cumSkippedWindows[i - 1];
        else
            cumSkippedWindows[i] = cumSkippedWindows[i - 1] + 1;
    }
    // Use the count matrix for resizing all new windows
    clear(c);
    resize(c, totalWinNum, Exact());
    unsigned skipped = 0;
    for (unsigned w = 0; w < totalWinNum; ++w)
    {
        c[w].chrom = o.chrom;
        while (!windowsWithContent[w + skipped])
            ++skipped;

        c[w].beginPos = beginPos  + (w + skipped) * newWinSize;
        resize(c[w].records, rgNum, Exact());
        for (unsigned rg = 0; rg < rgNum; ++rg)
        {
            const unsigned currentInsertCount = winCountMap[rg][w+skipped];
            if (currentInsertCount != 0)
            {
                reserve(c[w].records[rg], currentInsertCount, Exact());
            }
        }
    }
    // Now we can finally fill the new string of windows.
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        for (unsigned i = 0; i < length(o.records[rg]); ++i)
        {
            const unsigned currentWin = o.records[rg][i].pos;
            unsigned winIndex = (currentWin - beginPos) / newWinSize;
            // Adjust for the empty windows.
            SEQAN_ASSERT_GEQ(winIndex, 0u);
            SEQAN_ASSERT_GEQ(winIndex, cumSkippedWindows[winIndex]);
            winIndex -= cumSkippedWindows[winIndex];
            SEQAN_ASSERT_LT(winIndex, length(c));
            appendValue(c[winIndex].records[rg], o.records[rg][i]);
        }
    }
    return true;
}

#endif /* WINDOW_POPDEL_H_ */
