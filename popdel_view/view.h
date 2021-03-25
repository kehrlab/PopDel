#ifndef VIEW_POPDEL_H_
#define VIEW_POPDEL_H_

#include "popdel_view_parameter_parsing.h"
#include "../popdel_profile/window_podel.h"
#include "../utils_popdel.h"
#include "../insert_histogram_popdel.h"

using namespace seqan;

// =======================================================================================
// Function checkWindowInRegion()
// =======================================================================================
// Check if the given window is in the ROI
// Return -1 if the window is a on different contig or downstream of the ROI.
// Return 0 if the window is on the same chromosome but upstream of the ROI.
// Return 1 if the window is within the ROI.
inline int checkWindowInRegion(const Window & window,
                               const GenomicRegion & region,
                               const String<CharString> & contigNames)
{
    if (contigNames[window.chrom] != region.seqName ||
       (contigNames[window.chrom] == region.seqName &&
        window.beginPos >= region.endPos))
        return -1;

    else if (region.beginPos > window.beginPos)
        return 0;
    else
        return 1;
}
// =======================================================================================
// Function readAndWriteAllWindows()
// =======================================================================================
// Read and write all windows in the profile input stream.
template<typename TStream>
inline void readAndWriteAllWindows(TStream & in, const String<CharString> & contigNames, const unsigned & numReadGroups)
{
    Window window;
    while(readWindow(in, window, numReadGroups))
          writeWindow(std::cout, window, contigNames);
}
// Wrapper for handling the decompression of the input stream if necessary.
template<typename TStream>
inline void readAndWriteAllWindows(TStream & in,
                                   const String<CharString> & contigNames,
                                   const unsigned & numReadGroups,
                                   const bool & uncompressedIn)
{
    if (uncompressedIn)
    {
        readAndWriteAllWindows(in, contigNames, numReadGroups);
    }
    else
    {
        zlib_stream::zip_istream unzipper(in);
        readAndWriteAllWindows(unzipper, contigNames, numReadGroups);
    }
}

// =======================================================================================
// Function readAndWriteUncompressedWindows()
// =======================================================================================
// Reads the input windows and writes them without compression.
template<typename TStream>
inline void readAndWriteUncompressedWindows(TStream & in,
                                            std::ofstream & out,
                                            String<String<uint64_t> > & indexFields,
                                            const unsigned & indexRegionSize,
                                            const String<Histogram> & histograms)
{
    Window window;
    __int32 regionChrom = -1;
    __int32 regionPos = -1;
    zlib_stream::zip_istream unzipper(in);
    while(readWindow(unzipper, window, histograms))
    {
        if (indexNeedsUpdate(true, regionChrom, regionPos, window.chrom, window.beginPos, indexRegionSize))
            updateIndexFields(indexFields, regionChrom, regionPos, window.chrom, window.beginPos, out, indexRegionSize);

        writeWindow(out, window, histograms);
    }
}
// =======================================================================================
// Function readAndWriteUncompressedWindowsInRoi()
// =======================================================================================
// Reads the input windows and writes them without compression.
template<typename TStream>
inline void readAndWriteUncompressedWindowsInRoi(TStream & in,
                                                 std::ofstream & out,
                                                 String<String<uint64_t> > & indexFields,
                                                 const String<CharString> & contigNames,
                                                 const unsigned & indexRegionSize,
                                                 const String<Histogram> & histograms,
                                                 const GenomicRegion & region)
{
    Window window;
    __int32 regionChrom = -1;
    __int32 regionPos = -1;
    zlib_stream::zip_istream unzipper(in);
    while(readWindow(unzipper, window, histograms))
    {
        int retCode = checkWindowInRegion(window, region, contigNames);
        if (retCode == -1)              // Wrong contig or downstream of ROI
            break;
        else if (retCode == 0)          // Same contig, but upstream of ROI
            continue;
        else // (if retCode == 1)       // Window is in ROI     // TODO Continue here and check logic for updating index
        {
            if (indexNeedsUpdate(true, regionChrom, regionPos, window.chrom, window.beginPos, indexRegionSize))
                updateIndexFields(indexFields, regionChrom, regionPos, window.chrom, window.beginPos, out, indexRegionSize);

            writeWindow(out, window, histograms);
        }
    }
}
// =======================================================================================
// Function checkWindowInRegion()
// =======================================================================================
// Read the windows and print them if they are within the ROI
template<typename TStream>
inline void readAndWriteWindowsInRoi(TStream & in,
                                     const GenomicRegion & region,
                                     const String<CharString> & contigNames,
                                     const unsigned & numReadGroups)
{
    Window window;
    while(readWindow(in, window, numReadGroups))
    {
        int retCode = checkWindowInRegion(window, region, contigNames);
        if (retCode == -1)              // Wrong contig or downstream of ROI
            break;
        else if (retCode == 0)          // Same contig, but upstream of ROI
            continue;
        else // (if retCode == 1)       // Window is in ROI
        writeWindow(std::cout, window, contigNames);
    }
}
// Wrapper for handling decompression of the stream if necessary.
template<typename TStream>
inline void readAndWriteWindowsInRoi(TStream & in,
                                     const GenomicRegion & region,
                                     const String<CharString> & contigNames,
                                     const unsigned & numReadGroups,
                                     const bool & uncompressedIn)
{
    if (uncompressedIn)
    {
        readAndWriteWindowsInRoi(in, region, contigNames, numReadGroups);
    }
    else
    {
        zlib_stream::zip_istream unzipper(in);
        readAndWriteWindowsInRoi(unzipper, region, contigNames, numReadGroups);
    }
}
// =======================================================================================
// Function viewProfile()
// =======================================================================================
// Read the profile and print the content windows.
template <typename TStream>
inline int viewProfile(TStream & in,
                       const String<CharString> & contigNames,
                       const String<int32_t> & contigLengths,
                       const unsigned numReadGroups,
                       PopDelViewParameters & params)
{
    if (empty(params.regions))
    {
        readAndWriteAllWindows(in, contigNames, numReadGroups, params.uncompressedIn);
    }
    else
    {
        // Read file offset from index and move the stream there.
        fillInvalidPositions(params.regions[0], contigNames, contigLengths);
        jumpToRegion(in, contigNames, contigLengths, params.indexRegionSize, params.regions[0]);
        readAndWriteWindowsInRoi(in, params.regions[0], contigNames, numReadGroups, params.uncompressedIn);
    }
    return 0;
}
// =======================================================================================
// Function partitionChromosome()
// =======================================================================================
// Partition the given chromosome and store the regions
inline void partitionChromosome(PopDelViewParameters & params,
                                const CharString & contigName,
                                const unsigned & contigLength)
{
    unsigned numRegions = contigLength / params.chunkSize + 1;
    for (unsigned i = 0; i < numRegions; ++i)
    {
        GenomicRegion r;
        r.seqName = contigName;
        r.beginPos = i * params.chunkSize;
        r.endPos = std::min(i * params.chunkSize + params.chunkSize + params.padding, contigLength);
        appendValue(params.regions, r);
    }
}
// =======================================================================================
// Function partitionGenome()
// =======================================================================================
// Parition the whole genome and store the regions in params.regions.
// If a region is specified in params.regions[0] before, only that chromosome will be partiioned and the rest
//  of the profile will be ignored.
inline void partitionGenome(const String<CharString> & contigNames,
                            const String<unsigned> & contigLengths,
                            PopDelViewParameters & params)
{
    GenomicRegion region;
    if (!empty(params.regions))
        region = params.regions[0];

    clear(params.regions);
    unsigned totalRegionNum = 0;
    if (region.seqName == "")
    {
        for (unsigned i = 0; i < length(contigLengths); ++i)
            totalRegionNum += contigLengths[i] / params.chunkSize + 1;

        reserve(params.regions, totalRegionNum, Exact());
        for (unsigned i = 0; i < length(contigLengths); ++i)
            partitionChromosome(params, contigNames[i], contigLengths[i]);
    }
    else    // Limit partitioning to the given chromosome
    {
        unsigned rID = 0;
        for (; rID < length(contigLengths); ++rID)
        {
            if (region.seqName == contigNames[rID])
            {
                totalRegionNum += contigLengths[rID] / params.chunkSize + 1;
                break;
            }
        }
        reserve(params.regions, totalRegionNum, Exact());
        partitionChromosome(params, contigNames[rID], contigLengths[rID]);
    }
}
// =======================================================================================
// Function createChunkFileName()
// =======================================================================================
// Generate the filename based on the outfilename and the genomic region.
// Return the the name for the output file
inline void createChunkFileName(CharString & fileName,
                                const GenomicRegion & region,
                                const PopDelViewParameters & params)
{
        fileName = params.outfile;
        appendValue(fileName, '.');
        append(fileName, region.seqName);
        appendValue(fileName, '.');
        append(fileName, std::to_string(region.beginPos));
        appendValue(fileName, '-');
        append(fileName, std::to_string(region.endPos));
        append(fileName, ".profile");
}
// =======================================================================================
// Function unzipProfile()
// =======================================================================================
// Read the profile and save it uncompressed to the output file.
template <typename TStream>
inline int unzipProfile(TStream & in,
                        const String<CharString> & contigNames,
                        const String<int32_t> & contigLengths,
                        String<CharString> & readGroups,
                        String<Histogram> & histograms,
                        PopDelViewParameters & params)
{
    String<String<uint64_t> > indexFields;
    String<CharString> outFileNames;
    unsigned indexSize = resizeIndexFields(indexFields, contigLengths, params.outIndexRegionSize);
    if (params.chunkSize == 0) // No partioning desired.
    {
        std::ofstream out(toCString(params.outfile), std::ios::out | std::ios::binary);
        writeProfileHeader(out,
                          params.outIndexRegionSize,
                          indexSize,
                          readGroups,
                          histograms,
                          contigNames,
                          contigLengths,
                          false);
        if (empty(params.regions))  // Process the whole profile
        {
            readAndWriteUncompressedWindows(in,
                                            out,
                                            indexFields,
                                            params.outIndexRegionSize,
                                            histograms);
        }
        else    // Process a single region
        {
            fillInvalidPositions(params.regions, contigNames, contigLengths);
            jumpToRegion(in, contigNames, contigLengths, params.indexRegionSize, params.regions[0]);
            readAndWriteUncompressedWindowsInRoi(in,
                                                 out,
                                                 indexFields,
                                                 contigNames,
                                                 params.outIndexRegionSize,
                                                 histograms,
                                                 params.regions[0]);
        }
        writeIndexIntoHeader(out, indexFields, out.tellp());
    }
    else // Perform partitioning
    {
        partitionGenome(contigNames, contigLengths, params);
        std::ostringstream msg;
        msg << "Unzipping and partitioning the profile into " << length(params.regions) << " chunks of size " << params.chunkSize << "+" << params.padding << " bp: ";
        for (unsigned i = 0; i < length(params.regions); ++i)
        {
            msg << params.regions[i].seqName << ":" << params.regions[i].beginPos << "-" << params.regions[i].endPos;
            if (i + 1 != length(params.regions))
                msg << ", ";
        }
        printStatus(msg);
        for (unsigned r = 0; r < length(params.regions); ++r)
        {
            CharString outFileName;
            createChunkFileName(outFileName, params.regions[r], params);
            appendValue(outFileNames, outFileName);
            std::ofstream out(toCString(outFileName), std::ios::out | std::ios::binary);
            writeProfileHeader(out,
                               params.outIndexRegionSize,
                               indexSize,
                               readGroups,
                               histograms,
                               contigNames,
                               contigLengths,
                               false);
            jumpToRegion(in, contigNames, contigLengths, params.indexRegionSize, params.regions[r]);
            readAndWriteUncompressedWindowsInRoi(in,
                                                 out,
                                                 indexFields,
                                                 contigNames,
                                                 params.outIndexRegionSize,
                                                 histograms,
                                                 params.regions[r]);

            writeIndexIntoHeader(out, indexFields, out.tellp());
        }
    }
    std::ostringstream msg;
    if (params.chunkSize > 0)
    {
        msg << length(params.regions) << " unzipped profiles written to ";
        for (unsigned i = 0; i < length(outFileNames); ++i)
        {
            msg << "\'" << outFileNames[i] << "\'";
            if (i + 1 < length(outFileNames))
                msg << ", ";
        }
    }
    else
    {
        msg << "Unzipped profile written to \'" << params.outfile << "\'.";
    }
    printStatus(msg);
    return 0;
}

#endif /*VIEW_POPDEL_H_*/