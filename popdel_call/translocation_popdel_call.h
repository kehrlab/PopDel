#ifndef TRANSLOCATION_POPDEL_CALL_H_
#define TRANSLOCATION_POPDEL_CALL_H_

#include "../popdel_profile/profile_translocation_popdel.h"
#include "../popdel_profile/window_popdel.h"

// =======================================================================================
// Struct FirstTranslocationRead
// =======================================================================================
// Similar to TranclocationRead, but without member for refID because this info will be stored in the window.
// Should be used in the class TranslocationWindowEntry.
struct FirstTranslocationRead
{
    uint32_t        pos;            // Position on the reference
    unsigned char   clip_0;           // Number of soft clipped bases at the 5'-end of the read.
    unsigned char   clip_1;           // Number of soft clipped bases at the 3'-end of the read.

    FirstTranslocationRead(): pos(0), clip_0(0), clip_1(0) {}
};
// =======================================================================================
// Struct FirstTranslocationRead
// =======================================================================================
// Similar to TranslocationPair, but with a FirstTranslocationRead object for first to avoid storing of 
// redundant refID information in every first reads.
// Should be used in the class TranslocationWindow.
struct TranslocationWindowEntry
{
    FirstTranslocationRead first;
    TranslocationRead second;
    Orientation orientation;            // Orientation of first and second, in this order.

    TranslocationWindowEntry():
        first(), second(), orientation(Orientation::FR)
    {}

    TranslocationWindowEntry(unsigned p, unsigned c, Orientation o)
    {
        first.pos = p;
        first.clip_1 = c;
        orientation = o;
    }
    TranslocationWindowEntry(unsigned p, unsigned c_0, unsigned c_1, Orientation o)
    {
        first.pos = p;
        first.clip_0 = c_0;
        first.clip_1 = c_1;
        orientation = o;
    }
    TranslocationWindowEntry(const FirstTranslocationRead & f, const TranslocationRead & s, Orientation o)
    {
        first = f;
        second = s;
        orientation = o;
    }
};
// =======================================================================================
// Struct TranslocationWindow
// =======================================================================================
// Defines a translocationWindow window on the Genome. Can contain samples from multiple read groups.
struct TranslocationWindow
{
    int32_t chrom;
    int32_t beginPos;
    String<String<TranslocationWindowEntry> > records;

    TranslocationWindow():
        chrom(GenomicRegion::INVALID_ID), beginPos(GenomicRegion::INVALID_POS)
    {}

    TranslocationWindow(int32_t c, int32_t pos):
        chrom(c), beginPos(pos)
    {}
};
// =======================================================================================
// Struct TranslocationProfile
// =======================================================================================
// Class for buffering the windows of translocated records
struct TranslocationProfile
{
    int32_t chrom;
    unsigned currentPos;
    String<String<TranslocationWindowEntry> >  records;         // 1st index: Read group; 2nd index record
    String<Pair<Iterator<const String<TranslocationWindowEntry> >::Type> > activeReads; // On iterator pair per
                                                                                        // read group. (start, end);

    TranslocationProfile(const unsigned numReadGroups)
    {
        chrom = -1;
        currentPos = maxValue<unsigned>();
        resize(records, numReadGroups, Exact());
        resize(activeReads, numReadGroups, Exact());
    }
    // =======================================================================================
    // Struct add()
    // =======================================================================================
    // Add a record to the translocation profile for the given read group.
    inline void add(const unsigned readGroup, const TranslocationWindowEntry & e)
    {
        SEQAN_ASSERT_LT(readGroup, length(records));
        appendValue(records[readGroup], e);
    }
    // =======================================================================================
    // Struct getTotalTranslocationReadCounts()
    // =======================================================================================
    // Return the number of translocation reads in the current window
    // for all given RGs w/o considering contigs or orientations.
    inline unsigned getTotalTranslocationReadCounts(const TReadGroupIndices & rgs) const
    {
        unsigned c = 0;
        for(auto rg = begin(rgs); rg != end(rgs); ++rg)
            c += std::distance(activeReads[*rg].i1, activeReads[*rg].i2);

        // std::cout << "Found " << c << " translocation reads in window " << chrom << ":" << currentPos << std::endl;
        // std::cout << "Iterators are " << activeReads[0].i1  << " and " << activeReads[0].i2 << std::endl;
        return c;
    }
    // =======================================================================================
    // Struct clear()
    // =======================================================================================
    // Clear all records of all read groups in the translocation profile.
    inline void clear()
    {
        chrom = -1;
        currentPos = maxValue<unsigned>();
        for (unsigned rg = 0; rg < length(records); ++rg)
        {
            seqan::clear(records[rg]);
            activeReads[rg].i1 = end(records[rg]);
            activeReads[rg].i2 = end(records[rg]);
        }
    }
    // =======================================================================================
    // Struct setWindow()
    // =======================================================================================
    // Update the iterators in activeReads s.t. all pairs function as borders for the reads active in the given window.
    // If chrom is not equal to rID the check is performed starting from the start of records,
    // and from activeReads[rg] on otherwise.
    // Note: This function works best if the window is only incrementent in small steps because it is not
    // performing a binary search.
    inline void  setWindow(const int rID, const unsigned winStartPos, const unsigned windowSize = 30)
    {
        //std::cout << "CurrentPos is " << rID << ":" << winStartPos << "." << std::endl;
        if (winStartPos == 150)
            std::cout  << "";
       // std::cout << "Size of records:" <<length(records) << std::endl;
        for (unsigned rg = 0; rg < length(records); ++rg)
        {
            Iterator<const String<TranslocationWindowEntry> >::Type rec;
            if (chrom != rID)
            {
                //std::cout << rID << " != " << chrom << ". Starting from beginning of translocationRecords." << std::endl;
                rec = begin(records[rg]);
                chrom = rID;
            }
            else
            {
                rec = activeReads[rg].i1;       // TODO: Continue here and check why its pointing to 0x0 for rg=1
            }
            // Iterate and find first record after start of window
            while (rec != end(records[rg]) && rec->first.pos < winStartPos)
                ++rec;

            activeReads[rg].i1 = rec;
            // Iterate and find first record after end of window
            while (rec != end(records[rg]) && rec->first.pos < winStartPos + windowSize)
                ++rec;

            activeReads[rg].i2 = rec;
        }
        // For testing:
        if (activeReads[0].i1 != end(records[0]))
        {
            //std::cout << "Setting transloc window to " << chrom << ":" << activeReads[0].i1->first.pos << std::endl;
            //std::cout << "Iterators at " << activeReads[0].i1 << " and " << activeReads[0].i2 << std::endl;
        }
        else
        {
            //std::cout << "Records at " << chrom << " at end." << std::endl;
        }
        currentPos = winStartPos;
    }
};
// =======================================================================================
// Function readTranslocationCoordinates()
// =======================================================================================
// Read the chromosome and position of the translocation window from the input stream.
// Return true on success, false if EOF has been reached.
inline bool readTranslocationCoordinates(int32_t & chrom, int32_t & beginPos, zlib_stream::zip_istream & stream)
{
    // Read chromosome and position.
    stream.read(reinterpret_cast<char *>(&chrom), sizeof(uint32_t));
    if (stream.eof())
        return false;
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig name in translocation block."));
    stream.read(reinterpret_cast<char *>(&beginPos), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read position in translocation block."));
    return true;
}
// =======================================================================================
// Function unpackOrientations()
// =======================================================================================
// Put the orientations in groups of 4, s.t. they can fill all 8 bits of a uint8_t. Put these uints8_t into a String.
// Overload of unpackOrientations in window_popdel.h for translocation pairs. Does not rely on IDs.
inline void unpackOrientations(String<TranslocationWindowEntry> & records, const String<uint8_t> & packed)
{
    uint8_t mask = 192;     // Translates to 11000000 in binary format.
    for (unsigned i = 0; i < length(records); ++i)
        records[i].orientation = static_cast<Orientation>((packed[i / 4] & (mask >> i % 4 * 2)) >> (6-i%4*2));
}
// =======================================================================================
// Function readOrientations()
// =======================================================================================
// Read the packed string of read pair orientations from the stream and unpack them into the records.
// Overload of readOrientations in window.h for translocation windows
inline void readOrientations(String<TranslocationWindowEntry> & records,
                             zlib_stream::zip_istream & stream,
                             const unsigned numRecords)
{
    if (numRecords == 0)
        return;

    String<uint8_t> orientations;
    unsigned orientationChars = std::ceil(static_cast<double>(numRecords) / 4);
    resize(orientations, orientationChars, Generous()); // TODO: Why Generous()? Try Exact()!
    stream.read(reinterpret_cast<char *>(&orientations[0]), orientationChars);
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read orientations of read pairs in translocation block."));
    unpackOrientations(records, orientations);   // Unpack the orientations and store them in the window.
}
// =======================================================================================
// Function readRecords()
// =======================================================================================
// Read all records of one read group in the current window of the input stream. Used by readWindow().
inline void readRecords(String<TranslocationWindowEntry> & records,
                        zlib_stream::zip_istream & stream,
                        const unsigned windowBeginPos,
                        const unsigned numRecords)
{
    for (uint32_t i = 0; i < numRecords; ++i)
    {
        stream.read(reinterpret_cast<char *>(&records[i].first.pos), sizeof(char));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read position offset in translocation block."));
        records[i].first.pos += windowBeginPos;
        stream.read(reinterpret_cast<char *>(&records[i].first.clip_0), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read clipping of first translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].first.clip_1), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read clipping of first translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.refID), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read refID of second translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.pos), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read position of second translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.clip_0), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read clipping of second translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.clip_1), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read clipping of second translocation read."));
    }
}
// =======================================================================================
// Function readTranslocationWindow()
// =======================================================================================
// Reset the window and get the translocation records from the stream
// Return true on success, false if EOF is encountered
inline bool readTranslocationWindow(TranslocationWindow & window,
                                    zlib_stream::zip_istream & stream,
                                    const unsigned numReadGroups)
{
    reset(window);
    if (!readTranslocationCoordinates(window.chrom, window.beginPos, stream))
        return false;
    resize(window.records, numReadGroups);
    for (unsigned rg = 0; rg < numReadGroups; ++rg)
    {
        unsigned numRecords = readNumRecords(stream);
        resize(window.records[rg], numRecords);
        readOrientations(window.records[rg], stream, numRecords);
        readRecords(window.records[rg], stream, window.beginPos, numRecords);
    }
    return true;
}
// =======================================================================================
// Function translocationIndexBeginPos()
// =======================================================================================
// Return the file position where the translocation index starts.
// FilePos is calculated by taking the begin pos of the normal index and adding the 
// size of the normal index.
inline unsigned translocationIndexBeginPos(const unsigned numNormalIndexRegions)
{
    return indexBeginPos() + numNormalIndexRegions * sizeof(uint64_t);
}
// =======================================================================================
// Function jumpToTranslocRegion()
// =======================================================================================
// Jump to the first translocation record on chrom.
template<typename TStream>
inline void jumpToTranslocRegion(TStream & stream,
                                 const String<CharString> & contigNames,
                                 const CharString & chrom,
                                 const unsigned indexOffset)
{
    SEQAN_ASSERT(!stream.eof());
    // Find position of chrom in contigNames
    uint64_t i = 0;
    for (; i < length(contigNames) && chrom != contigNames[i]; ++i);

    // Calculate the position of chrom in the index and seek it
    uint64_t offset = 0;
    stream.seekg(indexOffset + i * sizeof(uint64_t));

    // Read the offset and move the stream.
    stream.read(reinterpret_cast<char *>(&offset), sizeof(uint64_t));
    stream.seekg(offset);
}
// =======================================================================================
// Function jumpToTranslocBlock()
// =======================================================================================
// Jump to the start of the translocation block (i.e. the first byte after the translocation guard)
template<typename TStream>
inline void jumpToTranslocBlock(TStream & stream, const unsigned indexOffset)
{
    // Move to the position of the first entry in the translocation index.
    uint64_t offset = 0;
    stream.seekg(indexOffset);

    // Read the offset and move the stream.
    stream.read(reinterpret_cast<char *>(&offset), sizeof(uint64_t));
    stream.seekg(offset);
}
// =======================================================================================
// Function readTillRoi()
// =======================================================================================
// Read all translocation entries until pos is reached, ignoring everything befor the ROI.
// Return 0 on success, 2 if a new chromsome has been reached before reaching pos and 3 at EOF.
inline unsigned readTillRoi(TranslocationWindow & window,
                            zlib_stream::zip_istream & unzipper,
                            const CharString & chrom,
                            const unsigned &  pos,
                            const String<CharString> & contigNames,
                            const TReadGroupIndices & rg)
{
    do
    {
        if (!readTranslocationWindow(window, unzipper, length(rg)))
        {
            return 3;               // EOF
        }
        if (contigNames[window.chrom] != chrom)
        {
            return 2;               // New chromosome
        }
    }
    while (static_cast<unsigned>(window.beginPos + 255) < pos);
    return 0;
}
// =======================================================================================
// Function convertTranslocationWindow()
// =======================================================================================
// Convert a translocation window of size "oldWinSize" to a string of windows with "newWinSize".
// Return true on success, false if there are no entries in the window.
inline bool convertTranslocationWindow(String<TranslocationWindow> & c,
                                       const TranslocationWindow & o,
                                       unsigned oldWinSize = 256,
                                       unsigned newWinSize = 30)
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
            if (o.records[rg][0].first.pos < beginPos)
                beginPos = o.records[rg][0].first.pos ;
        }
    }
    beginPos = (beginPos / newWinSize) * newWinSize;

    String<bool> windowsWithContent;
    resize(windowsWithContent, std::ceil(static_cast<double>(oldWinSize + newWinSize) / newWinSize), false, Exact());
    // Now count how many entries there will be in each window for each read group.
    // Also keep track of the number of different windows.
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        for (unsigned i = 0; i < length(o.records[rg]); ++i)
        {
            const unsigned currentWin = o.records[rg][i].first.pos;
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
    // Create a string indicating how many windows were skipped before a given window.
    // E.g. winCountMap hold 9 windows, but 4,6,7 (starting from 0) are empty.
    // Thus cumWindowsWithContent will hold [0, 0, 0, 0, 0, 1, 1 , 2, 3]
    // We need this for adjusting the index when saving the entries in c.
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
            const unsigned currentWin = o.records[rg][i].first.pos;
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
// =======================================================================================
// Function addRgTranslocationRecordsToProfile()
// =======================================================================================
// Add all records of the given read groups to the TranslocationProfile.
inline void addRgTranslocationRecordsToProfile(TranslocationProfile & profile,
                                               const TReadGroupIndices & rgs,
                                               const TranslocationWindow & win)
{
    for (unsigned r = 0; r < length(win.records); ++r)
    {   // Add all records of this read group to the profile.
        for (unsigned i = 0; i < length(win.records[r]); ++i)
        {
            TEntryIdType rg = rgs[r];
            profile.add(rg, win.records[r][i]);
        }
    }
}
// =======================================================================================
// Function readTranslocationSegment()
// =======================================================================================
// Load the translocation records for the given ROI into the TranslocationProfile
inline unsigned readTranslocationSegment(TranslocationProfile & translocProfile,
                                         std::ifstream & file,
                                         const CharString & fileName,
                                         const TReadGroupIndices & rg,
                                         const String<CharString> & contigNames,
                                         const GenomicRegion & roi)
{
    zlib_stream::zip_istream unzipper(file);
    TranslocationWindow window;
    String<TranslocationWindow> convertedWindows;
    translocProfile.clear();
    while (true)
    {
        unsigned ret = readTillRoi(window, unzipper, roi.seqName, roi.beginPos, contigNames,  rg);
        if (ret > 0)
            return ret;

        if(!convertTranslocationWindow(convertedWindows, window))
        {
            std::cerr << "[PopDel] Error in profile \"" << fileName << "\"." << std::endl;
            std::ostringstream msg;
            msg << "[PopDel] Could not convert translocation window at " << contigNames[window.chrom] <<
            ":" << window.beginPos << ". The profile \"" << fileName << "\" might be corrupted.";
            SEQAN_THROW(IOError(toCString(msg.str())));
        }

        // for (unsigned i = 0; i < length(convertedWindows); ++i)  // Output for testing
        // {
        //     std::cout << convertedWindows[i].chrom << ":" << convertedWindows[i].beginPos << "\t";
        //     for(unsigned j = 0; j < length(convertedWindows[i].records); ++j)
        //         std::cout << convertedWindows[i].records[0][j].first.pos << ";";    // Only checking for RG 0

        //     std::cout << std::endl;
        // }
        for (auto it = begin(convertedWindows); it != end(convertedWindows); ++it)
            addRgTranslocationRecordsToProfile(translocProfile, rg,  *it);
    }
}
#endif /* TRANSLOCATION_POPDEL_CALL_H_ */