#ifndef WINDOW_POPDEL_H_
#define WINDOW_POPDEL_H_

#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include "../insert_histogram_popdel.h"

using namespace seqan;

// =======================================================================================
// Struct Window
// =======================================================================================
// Defines a window on the Genome. Can contain samples from multiple read groups.

struct Window
{
    typedef Pair<unsigned, unsigned> TRecord; // the end position of the forward read and the insert size.
    __int32 chrom;
    __int32 beginPos;
    String<String<TRecord> > insertSizes;

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
// Return true if the window is empty (i.e window.insertSizes has length 0)

inline bool isEmpty(const Window & window)
{
    return length(window.insertSizes) == 0;
}

// =======================================================================================
// Function reset()
// =======================================================================================
// Reset chrom, beginPos and insertSizes window but keep the same number of read groups.

inline void reset(Window & window)
{
    window.chrom = 0;
    window.beginPos = 0;
    for (unsigned rg = 0; rg < length(window.insertSizes); ++rg)
        clear(window.insertSizes[rg]);
}

// =======================================================================================
// Function addRecord()
// =======================================================================================
// Aggregate insertSize and numWindows as pair and add it to window for read group.

inline void addRecord(Window & window,
                      __int32 beginPos,
                      unsigned insertSize,
                      unsigned readGroup,
                      unsigned numReadGroups)
{
    SEQAN_ASSERT_NEQ(numReadGroups, 0u);
    SEQAN_ASSERT_LT(readGroup, numReadGroups);
    if (isEmpty(window))                                // TODO: Maybe move this to some initialization step.
        resize(window.insertSizes, numReadGroups);
    appendValue(window.insertSizes[readGroup], Pair<unsigned>(beginPos, insertSize));
    SEQAN_ASSERT_EQ(isEmpty(window), false);
}

// =======================================================================================
// Function writeWindow()
// =======================================================================================
// Write info from a window to output stream. Write CHROM, POS and insert size deviation from median per read group.

template<typename TStream>
inline void writeWindow(TStream & stream,
                        const Window & window,
                        const String<CharString> & contigNames,
                        const String<Histogram> & histograms)
{
    typedef typename Iterator<const String<Window::TRecord> >::Type TValueIter;
    SEQAN_ASSERT(!isEmpty(window));

    // Write chromosome and position.
    stream << contigNames[window.chrom];
    stream << "\t" << window.beginPos + 1;

    for (unsigned i = 0; i < length(histograms); ++i)
    {
        TValueIter it = begin(window.insertSizes[i]);
        TValueIter itEnd = end(window.insertSizes[i]);
        if (it == itEnd)
            stream << "\t" << ".";                                                  // No insert size in list.
        else                                                                    // Write the first insert size.
        {
            stream << "\t" << static_cast<int>((*it).i1 - window.beginPos);
            stream << ":" << static_cast<int16_t>((*it).i2 - _round(histograms[i].median));
            ++it;
        }

        while (it != itEnd)                                                     // Write the remaining insert sizes.
        {
            stream << "," << static_cast<int>((*it).i1 - window.beginPos);
            stream << ":" << static_cast<int16_t>((*it).i2 - _round(histograms[i].median));
            ++it;
        }
    }
    stream << std::endl;
}

template<typename TStream>
inline void writeWindow(TStream & stream,
                        const Window & window,
                        const String<CharString> & contigNames)
{
    typedef typename Iterator<const String<Window::TRecord> >::Type TValueIter;
    SEQAN_ASSERT(!isEmpty(window));

    // Write chromosome and position.
    stream << contigNames[window.chrom];
    stream << "\t" << window.beginPos + 1;

    for (unsigned i = 0; i < length(window.insertSizes); ++i)
    {
        TValueIter it = begin(window.insertSizes[i]);
        TValueIter itEnd = end(window.insertSizes[i]);
        if (it == itEnd)
            stream << "\t" << ".";                                                  // No insert size in list.
        else                                                                    // Write the first insert size.
        {
            stream << "\t" << static_cast<int>((*it).i1 - window.beginPos);
            stream << ":" << static_cast<int16_t>((*it).i2);
            ++it;
        }

        while (it != itEnd)                                                     // Write the remaining insert sizes.
        {
            stream << "," << static_cast<int>((*it).i1 - window.beginPos);
            stream << ":" << static_cast<int16_t>((*it).i2);
            ++it;
        }
    }
    stream << std::endl;
}

inline void writeWindow(zlib_stream::zip_ostream & stream,
                        const Window & window,
                        const String<Histogram> & histograms)
{
    typedef Iterator<const String<Window::TRecord> >::Type TValueIter;
    SEQAN_ASSERT(!isEmpty(window));

    // Write chromosome and position.
    uint32_t chrom = window.chrom;
    stream.write(reinterpret_cast<char *>(&chrom), sizeof(uint32_t));
    uint32_t beginPos = window.beginPos;
    stream.write(reinterpret_cast<char *>(&beginPos), sizeof(uint32_t));

    SEQAN_ASSERT_EQ(length(window.insertSizes), length(histograms));

    for (unsigned rg = 0; rg < length(window.insertSizes); ++rg)
    {
        // Write number of read pairs for this read group.
        uint32_t numInsertSizes = length(window.insertSizes[rg]);
        stream.write(reinterpret_cast<char *>(&numInsertSizes), sizeof(uint32_t));

        TValueIter it = begin(window.insertSizes[rg]);
        TValueIter itEnd = end(window.insertSizes[rg]);

        while (it != itEnd)                                                         // Write the read pairs.
        {
            if (static_cast<unsigned>((*it).i1) > beginPos + 255u ||
                static_cast<unsigned>((*it).i1) < beginPos)
                SEQAN_THROW(ParseError("Too big offset in Window!"));
            char posOffset = (*it).i1 - window.beginPos;
            stream.write(reinterpret_cast<char *>(&posOffset), sizeof(char));
            int16_t insertSizeDeviation = (*it).i2 - _round(histograms[rg].median);
            stream.write(reinterpret_cast<char *>(&insertSizeDeviation), sizeof(int16_t));
            ++it;
        }
    }
}
// Return the stating position of the window that contains the given position.
inline unsigned posToWin(const unsigned & pos, const unsigned windowSize = 30u)
{
    return (pos / windowSize) * windowSize;
}
// =======================================================================================
// Function readWindow()
// =======================================================================================
// Read the entires of the next window in theopened and unzipped stream.
// Note that this will potentially safe an signed value in an unsigned variable (window.i2).
// The correct value can be retrieved by casting this value to int16_t.
// Return true on sucess and false if EOF has been reached.
inline bool readWindow(zlib_stream::zip_istream & stream,
                       Window & window,
                       unsigned numReadGroups)
{
    reset(window);

    // Read chromosome and position.
    stream.read(reinterpret_cast<char *>(&window.chrom), sizeof(uint32_t));
    if (stream.eof())
        return false;
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig name."));
    stream.read(reinterpret_cast<char *>(&window.beginPos), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig position."));

    resize(window.insertSizes, numReadGroups);

    for (unsigned rg = 0; rg < numReadGroups; ++rg)
    {
        // Read the number of read pairs for this read group.
        uint32_t numInsertSizes = 0;
        stream.read(reinterpret_cast<char *>(&numInsertSizes), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read number of read pairs."));
        resize(window.insertSizes[rg], numInsertSizes);

        // Read the read pair data for this read group.
        for (uint32_t i = 0; i < numInsertSizes; ++i)
        {
            stream.read(reinterpret_cast<char *>(&window.insertSizes[rg][i].i1), sizeof(char));
            if (!stream.good())
                SEQAN_THROW(ParseError("[PopDel] Unable to read position offset."));
            window.insertSizes[rg][i].i1 += window.beginPos;
            stream.read(reinterpret_cast<char *>(&window.insertSizes[rg][i].i2), sizeof(int16_t));
            if (!stream.good())
                SEQAN_THROW(ParseError("[PopDel] Unable to read insert size deviation."));
        }
    }
    return true;
}
// Overload that adds the median insert-size to the deviation before storing it in window.
inline bool readWindow(zlib_stream::zip_istream & stream,
                       Window & window,
                       const String<Histogram> & histograms)
{
    reset(window);

    // Read chromosome and position.
    stream.read(reinterpret_cast<char *>(&window.chrom), sizeof(uint32_t));
    if (stream.eof())
        return false;
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig name."));
    stream.read(reinterpret_cast<char *>(&window.beginPos), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(ParseError("[PopDel] Unable to read contig position."));

    resize(window.insertSizes, length(histograms));

    for (unsigned rg = 0; rg < length(histograms); ++rg)
    {
        // Read the number of read pairs for this read group.
        uint32_t numInsertSizes = 0;
        stream.read(reinterpret_cast<char *>(&numInsertSizes), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(ParseError("[PopDel] Unable to read number of read pairs."));
        resize(window.insertSizes[rg], numInsertSizes);

        // Read the read pair data for this read group.
        for (uint32_t i = 0; i < numInsertSizes; ++i)
        {
            stream.read(reinterpret_cast<char *>(&window.insertSizes[rg][i].i1), sizeof(char));
            if (!stream.good())
                SEQAN_THROW(ParseError("[PopDel] Unable to read position offset."));
            window.insertSizes[rg][i].i1 += window.beginPos;
            stream.read(reinterpret_cast<char *>(&window.insertSizes[rg][i].i2), sizeof(int16_t));
            if (!stream.good())
                SEQAN_THROW(ParseError("[PopDel] Unable to read insert size deviation."));
            window.insertSizes[rg][i].i2 += _round(histograms[rg].median);
        }
    }
    return true;
}
// =======================================================================================
// Function convertWindow()
// =======================================================================================
// Convert a window of size "oldWinSize" to a string of windows with "newWinSize".  //TODO: Adapt to sorted orig. win.
// Return true on success, false if there are no entries in the window.
inline bool convertWindow(Window & o, String<Window> & c, unsigned oldWinSize = 256, unsigned newWinSize = 30)
{
    SEQAN_ASSERT_LT(newWinSize, oldWinSize);
    const unsigned rgNum = length(o.insertSizes);
    // Maxtrix used for counting how many elements will be in each window. First index: RG, second index: Window
    String<String<unsigned> > winCountMap;
    resize(winCountMap, rgNum, Exact());
    for (unsigned rg = 0; rg < rgNum; ++rg)
        resize(winCountMap[rg], std::ceil(static_cast<double>(oldWinSize + newWinSize) / newWinSize ), 0 , Exact());

    // Get the lowest starting position of all read groups
    unsigned beginPos = maxValue<unsigned>();
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        if (length(o.insertSizes[rg]) != 0)
        {
            if (o.insertSizes[rg][0].i1 < beginPos)
                beginPos = o.insertSizes[rg][0].i1 ;
        }
    }
    beginPos = (beginPos / newWinSize) * newWinSize;

    String<bool> windowsWithContent;
    resize(windowsWithContent, std::ceil(static_cast<double>(oldWinSize + newWinSize) / newWinSize),false, Exact());
    // Now count how many entries there will be in each window for each read group.
    // Also keep track of the number of different windows.
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        for (unsigned i = 0; i < length(o.insertSizes[rg]); ++i)
        {
            const unsigned currentWin = o.insertSizes[rg][i].i1;
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
        resize(c[w].insertSizes, rgNum, Exact());
        for (unsigned rg = 0; rg < rgNum; ++rg)
        {
            const unsigned currentInsertCount = winCountMap[rg][w+skipped];
            if (currentInsertCount != 0)
            {
                reserve(c[w].insertSizes[rg], currentInsertCount, Exact());
            }
        }
    }
    // Now we can finally fill the new string of windows.
    for (unsigned rg = 0; rg < rgNum; ++rg)
    {
        for (unsigned i = 0; i < length(o.insertSizes[rg]); ++i)
        {
            const unsigned currentWin = o.insertSizes[rg][i].i1;
            unsigned winIndex = (currentWin - beginPos) / newWinSize;
            // Adjust for the empty windows.
            SEQAN_ASSERT_GEQ(winIndex, 0u);
            SEQAN_ASSERT_GEQ(winIndex, cumSkippedWindows[winIndex]);
            winIndex -= cumSkippedWindows[winIndex];
            SEQAN_ASSERT_LT(winIndex, length(c));
            appendValue(c[winIndex].insertSizes[rg], o.insertSizes[rg][i]);
        }
    }
    return true;
}

#endif /* WINDOW_POPDEL_H_ */
