#ifndef POPDEL_VIEW_POPDEL_H_
#define POPDEL_VIEW_POPDEL_H_

#include "../popdel_profile/window_popdel.h"
#include "popdel_view_parameter_parsing.h"
#include "../popdel_call/translocation_popdel_call.h"


// =======================================================================================
// Function printWindow()
// =======================================================================================
// Print info from a window to the output stream. 
template<typename TStream>
inline void printWindow(TStream & stream,
                        const Window & window,
                        const String<CharString> & contigNames,
                        const bool printClipping,
                        const bool printOrientation)
{
    typedef typename Iterator<const String<ReadPair> >::Type TValueIter;
    SEQAN_ASSERT(!isEmpty(window));

    stream << contigNames[window.chrom];                                        // Write chromosome and position.
    stream << "\t" << window.beginPos + 1;

    for (unsigned i = 0; i < length(window.records); ++i)
    {
        TValueIter it = begin(window.records[i]);
        TValueIter itEnd = end(window.records[i]);
        if (it == itEnd)
            stream << "\t" << ".";                                              // No tip distance in list.
        else                                                                    // Write the first tip distance
        {
            stream << "\t" << static_cast<int>(it->pos - window.beginPos);
            stream << ":" << it->distance;
            if (printClipping)
                stream << "(" << it->clipping << ")";
            if (printOrientation)
                stream << ":" << it->orientation;
            ++it;
        }
        while (it != itEnd)                                                     // Write the remaining tip distances
        {
            stream << "," << static_cast<int>(it->pos - window.beginPos);
            stream << ":" << it->distance;
            if (printClipping)
                stream << "(" << it->clipping << ")";
            if (printOrientation)
                stream << ":" << it->orientation;
            ++it;
        }
    }
    stream << std::endl;
}
// =======================================================================================
// Function printWindows()
// =======================================================================================
// Decompress the stream and print all windows
template<typename TStream>
inline void printWindows(TStream & in,
                         const PopDelViewParameters & params,
                         const String<CharString> & contigNames,
                         const unsigned numReadGroups)
{
    Window window;
    zlib_stream::zip_istream unzipper(in);
    while(readWindow(window, unzipper, numReadGroups))
          printWindow(std::cout, window, contigNames, params.printClipping, params.printOrientation);
}
// =======================================================================================
// Function printWindowsInRegion()
// =======================================================================================
// Decompress the stream and print the windows in the region.
template<typename TStream>
inline void printWindowsInRegion(TStream & in,
                                 const PopDelViewParameters & params,
                                 const String<CharString> & contigNames,
                                 const unsigned numReadGroups)
{
    Window window;
    zlib_stream::zip_istream unzipper(in);
    while(readWindow(window, unzipper, numReadGroups))
    {
        if (contigNames[window.chrom] != params.region.seqName ||
            (contigNames[window.chrom] == params.region.seqName && window.beginPos >= params.region.endPos))
            break;

        if (params.region.beginPos > window.beginPos)
            continue;

        printWindow(std::cout, window, contigNames, params.printClipping, params.printOrientation);
    }
}

// =======================================================================================
// Function printTranslocationWindow()
// =======================================================================================
// Write info from a translocation window to the output stream.
// Note: Requires the stream to be already decompressed.
template<typename TStream>
inline void printTranslocationWindow(TStream & stream,
                                     const TranslocationWindow & window,
                                     const String<CharString> & contigNames,
                                     const bool printClipping,
                                     const bool printOrientation)
{
    typedef typename Iterator<const String<TranslocationWindowEntry> >::Type TValueIter;
    SEQAN_ASSERT(!isEmpty(window));

    stream << contigNames[window.chrom];                                        // Write chromosome and position.
    stream << "\t" << window.beginPos + 1;

    for (unsigned i = 0; i < length(window.records); ++i)
    {
        TValueIter itBegin = begin(window.records[i]);
        TValueIter it = itBegin;
        TValueIter itEnd = end(window.records[i]);
        if (it == itEnd)
            stream << "\t" << ".";                                              // Empty list.
        else
            stream << "\t";

        while (it != itEnd)                                                     // Write all entires
        {
            if (it != itBegin)
            stream << ",";
            stream << static_cast<int>(it->first.pos - window.beginPos);
            if (printClipping)
                stream << "(" << static_cast<unsigned>(it->first.clip) << ")";
            if (printOrientation)
            {
                stream << ":";
                if (firstReadForward(it->orientation))
                    stream << "F";
                else
                    stream << "R";
            }
            stream << "+" << contigNames[it->second.refID] << ":" << it->second.pos + 1;
            if (printClipping)
                stream << "(" << static_cast<unsigned>(it->second.clip) << ")";
            if (printOrientation)
            {
                stream << ":";
                if (secondReadForward(it->orientation))
                    stream << "F";
                else
                    stream << "R";
            }
            ++it;
        }
    }
    stream << std::endl;
}
// =======================================================================================
// Function printTranslocationWindows()
// =======================================================================================
// Decompress the stream and print all translocation windows
template<typename TStream>
inline void printTranslocationWindows(TStream & in,
                                      const PopDelViewParameters & params,
                                      const String<CharString> & contigNames,
                                      const unsigned numReadGroups)
{
    TranslocationWindow transWindow;
    zlib_stream::zip_istream unzipper(in);
    while(readTranslocationWindow(transWindow, unzipper, numReadGroups))
        printTranslocationWindow(std::cout,
                                 transWindow,
                                 contigNames,
                                 params.printClipping,
                                 params.printOrientation);
}
// =======================================================================================
// Function printTranslocationWindowsInRegion()
// =======================================================================================
// Decompress the stream and print the translocation windows in the given region.
template<typename TStream>
inline void printTranslocationWindowsInRegion(TStream & in,
                                              const PopDelViewParameters & params,
                                              const String<CharString> & contigNames,
                                              const unsigned numReadGroups)
{
    TranslocationWindow transWindow;
    zlib_stream::zip_istream unzipper(in);
    while(readTranslocationWindow(transWindow, unzipper, numReadGroups))
    {
        if (contigNames[transWindow.chrom] != params.region.seqName ||
            (contigNames[transWindow.chrom] == params.region.seqName && transWindow.beginPos >= params.region.endPos))
            break;
        if (params.region.beginPos > transWindow.beginPos)
            continue;

        printTranslocationWindow(std::cout,
                                 transWindow,
                                 contigNames,
                                 params.printClipping,
                                 params.printOrientation);
    }
}
#endif /* POPDEL_VIEW_PARAMETER_PARSING_POPDEL_H_*/