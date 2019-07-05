#ifndef INSERT_HISTOGRAM_POPDEL_H_
#define INSERT_HISTOGRAM_POPDEL_H_

#include <math.h>

#include <seqan/sequence.h>

#include "utils_popdel.h"

using namespace seqan;

// -----------------------------------------------------------------------------
// Typedefs
// -----------------------------------------------------------------------------
typedef Iterator<String<double>, Rooted>::Type THistIter;
// -----------------------------------------------------------------------------
// Struct Histogram
// -----------------------------------------------------------------------------
struct Histogram
{                          // Counts of insert sizes. Index 1 holds the counts of inserts of len. 1 (+offset) etc
    String<double> values; // First and last index are reserved for irregular insert sizes (unmapped pair, to big, etc.)
    double         min_prob;               // 1/10000 of the maximum count in "values"
    int            offset;                 // Difference of element index 0 in the histogram string from insert size 0
    double         mean;
    double         stddev;
    unsigned       median;                 // Index of median (~median insert size)
    unsigned       readLength;             // length of one read making up the read pair
    double         coverage;
    unsigned       lowerQuantileDist;      // Distance from median to 5% quantile.
    unsigned       upperQuantileDist;      // Distance from median to 95% quantile.
    unsigned       windowSize;             // The size of one window in bp's.

    Histogram() :
    values(""),
    min_prob(0.0),
    offset(0),
    mean(0.0),
    stddev(0.0),
    median(0),
    readLength(0),
    coverage(0),
    windowSize(0)
    {}
    // ---------------------------------------------------------------------------------------
    // Function clearStrings()
    // ---------------------------------------------------------------------------------------
    // Clears all member strings.
    inline void clearStrings()
    {
        clear(values);
    }
};
inline unsigned getHistLeftBorder(const Histogram & hist)
{
    SEQAN_ASSERT_GT(hist.stddev, 0);
    return std::max(1, static_cast<int>(std::floor(static_cast<int>(hist.median) - 3 * hist.stddev)));
}
inline unsigned getHistRightBorder(const Histogram & hist)
{
    SEQAN_ASSERT_GT(length(hist.values), 0u);
    SEQAN_ASSERT_GT(hist.stddev, 0);
    return std::min(length(hist.values) - 1, (size_t)std::ceil(hist.median + 3 * hist.stddev) + 1);
}
// =======================================================================================
// Function writeProfileHeader()
// =======================================================================================
// Write the header (CHROM POS RG1 RG2 ... RGn) to the output stream.

template<typename TStream>
inline void writeProfileHeader(TStream & stream,
                               const String<CharString> & readGroups,
                               const String<CharString> & contigNames,
                               const String<int32_t> & contigLengths)
{
    // Write version line.
    stream << "@HD";
    stream << "\t" << "VN:" << "0.1";
    stream << std::endl;

    // Write read group names.
    for (unsigned i = 0; i < length(readGroups); ++i)
    {
        stream << "@RG";
        stream << "\t" << "ID:" << readGroups[i];
        stream << std::endl;
    }

    // Write chromosome names.
    for (unsigned i = 0; i < length(contigNames); ++i)
    {
        stream << "@SQ";
        stream << "\t" << "SN:" << contigNames[i];
        stream << "\t" << "LN:" << contigLengths[i];
        stream << std::endl;
    }
}

template<typename TStream>
inline void writeProfileHeader(TStream & stream,
                               const std::map<CharString, unsigned> & readGroups,
                               const String<CharString> & contigNames,
                               const String<int32_t> & contigLengths)
{
    typedef typename std::map<CharString, unsigned>::const_iterator TIter;

    // Get the IDs of all read groups.
    String<CharString> rg;
    resize(rg, length(readGroups));
    for (TIter it = readGroups.begin(); it != readGroups.end(); ++it)
        rg[it->second] = it->first;

    writeProfileHeader(stream, rg, contigNames, contigLengths);
}

template<typename TStream>
inline void writeProfileHeader(TStream & stream,
                               const unsigned & regionSize,
                               const unsigned & numRegions,
                               const std::map<CharString, unsigned> & readGroups,
                               String<Histogram> & histograms,
                               const String<CharString> & contigNames,
                               const String<int32_t> & contigLengths,
                               const bool mergeRG)
{
    typedef typename std::map<CharString, unsigned>::const_iterator TIter;

    SEQAN_ASSERT_EQ(length(contigNames), length(contigLengths));

    // Get the IDs of all read groups.
    String<CharString> rg;
    if (mergeRG)
    {
        resize(rg, 1, Exact());
        rg[0] = readGroups.begin()->first;
    }
    else
    {
        resize(rg, length(readGroups));
        for (TIter it = readGroups.begin(); it != readGroups.end(); ++it)
            rg[it->second] = it->first;
    }

    // Write magic string.
    stream.write("POPDEL\1", 7);

    // Write index regions size and size of the index.
    stream.write(reinterpret_cast<const char*>(&regionSize), sizeof(uint32_t));
    stream.write(reinterpret_cast<const char *>(&numRegions), sizeof(uint32_t));

    // Write index spaceholder.
    uint64_t pos = 0;
    for (unsigned i = 0; i < numRegions; ++i)
        stream.write(reinterpret_cast<char *>(&pos), sizeof(uint64_t));

    // Write number of read groups.
    uint32_t numReadGroups = mergeRG?1:readGroups.size();
    stream.write(reinterpret_cast<char *>(&numReadGroups), sizeof(uint32_t));

    // Write read group names and insert size histograms.
    for (unsigned i = 0; i < numReadGroups; ++i)
    {
        uint32_t rgLen = length(rg[i]) + 1;
        stream.write(reinterpret_cast<char *>(&rgLen), sizeof(uint32_t));
        stream.write(reinterpret_cast<char *>(&rg[i][0]), length(rg[i]));
        stream.write("\0", 1);
        stream.write(reinterpret_cast<char *>(&histograms[i].median), sizeof(uint16_t));
        stream.write(reinterpret_cast<char *>(&histograms[i].stddev), sizeof(double));
        stream.write(reinterpret_cast<char *>(&histograms[i].readLength), sizeof(uint16_t));
        unsigned histStart = getHistLeftBorder(histograms[i]);
        unsigned histEnd = getHistRightBorder(histograms[i]);
        stream.write(reinterpret_cast<char *>(&histStart), sizeof(uint16_t));
        stream.write(reinterpret_cast<char *>(&histEnd), sizeof(uint16_t));
        THistIter it = begin(histograms[i].values) + histStart;
        THistIter itEnd = begin(histograms[i].values) + histEnd;
        for (; it != itEnd; ++it)
        {
            double val = *it;
            stream.write(reinterpret_cast<char *>(&val), sizeof(double));
        }
    }

    // Write number of chromosomes.
    uint32_t numContigs = length(contigNames);
    stream.write(reinterpret_cast<char *>(&numContigs), sizeof(uint32_t));

    // Write chromosome names and lengths.
    for (unsigned i = 0; i < numContigs; ++i)
    {
        uint32_t contigNameLen = length(contigNames[i]) + 1;
        stream.write(reinterpret_cast<char *>(&contigNameLen), sizeof(uint32_t));
        stream.write(reinterpret_cast<const char *>(&contigNames[i][0]), length(contigNames[i]));
        stream.write("\0", 1);
        stream.write(reinterpret_cast<const char *>(&contigLengths[i]), sizeof(uint32_t));
    }
}

// =======================================================================================
// Function writeIndexIntoHeader()
// =======================================================================================
template<typename TStream>
inline void writeIndexIntoHeader(TStream & stream, String<String<uint64_t> > profileIndex, uint64_t maxOffset)
{
    // Fill in file offsets for empty regions.
    uint64_t prev = maxOffset;
    for (int i = length(profileIndex) - 1; i >= 0; --i)
    {
        for (int j = length(profileIndex[i]) - 1; j >= 0; --j)
        {
            if (profileIndex[i][j] == 0)
                profileIndex[i][j] = prev;
            else
                prev = profileIndex[i][j];
        }
    }

    // Save current position of stream and move it to index begin pos.
    uint64_t last = stream.tellp();
    stream.seekp(7 + 2 * sizeof(uint32_t));

    // Write the index.
    for (unsigned i = 0; i < length(profileIndex); ++i)
    {
        for (unsigned j = 0; j < length(profileIndex[i]); ++j)
        {
            stream.write(reinterpret_cast<char *>(&profileIndex[i][j]), sizeof(uint64_t));
        }
    }

    // Move current position of stream back to ending.
    stream.seekp(last); 
}

// =======================================================================================
// Function readProfileHeader()
// =======================================================================================
template<typename TStream>
inline void readProfileHeader(TStream & stream,
                              const CharString & filename,
                              String<CharString> & readGroups,
                              String<Histogram> & histograms,
                              String<CharString> & contigNames,
                              String<int32_t> & contigLengths,
                              unsigned & indexRegionSize)
{
    // Read the magic string.
    CharString buffer;
    resize(buffer, 7);
    stream.read(&buffer[0], 7);

    //Prepare message for potential errors.
    std::ostringstream msg;
    msg << "[PopDel] Corrupted profile \"" << filename << "\": ";

    if (!stream.good())
    {
        msg << "Unable to read magic string.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
    if (buffer != "POPDEL\1")
    {
        msg << "Magic string is wrong.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }

    // Read the index regions size, size of the index and skip (ignore) the index.
    stream.read(reinterpret_cast<char *>(&indexRegionSize), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read index region size.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
    unsigned numRegions = 0;
    stream.read(reinterpret_cast<char *>(&numRegions), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read index size.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
    stream.ignore(numRegions * sizeof(uint64_t));

    // Read the number of read groups.
    uint32_t numReadGroups = 0;
    stream.read(reinterpret_cast<char *>(&numReadGroups), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read number of read groups.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
    if (numReadGroups == 0) // How can this even happen?
    {
        msg << "Profile seems to have 0 read groups.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
    resize(readGroups, numReadGroups);
    resize(histograms, numReadGroups);

    // Read the read group names and histograms.
    for (unsigned i = 0; i < numReadGroups; ++i)
    {
        // Read length of read group name.
        uint32_t nameLen = 0;
        stream.read(reinterpret_cast<char *>(&nameLen), sizeof(uint32_t));
        if (!stream.good())
        {
            msg << "Unable to read length of read group name.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }

        // Read the read group name.
        resize(readGroups[i], nameLen-1);
        stream.read(reinterpret_cast<char *>(&readGroups[i][0]), nameLen-1);
        stream.read(&buffer[0], 1);
        if (!stream.good())
        {
            msg << "Unable to read read group name.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        if (buffer[0] != '\0') // Expect read group names to be null-terminated.
        {
            msg << "Expecting read group name to be null-terminated.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        // Read the insert size statistics.
        stream.read(reinterpret_cast<char *>(&histograms[i].median), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read median insert size.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        stream.read(reinterpret_cast<char *>(&histograms[i].stddev), sizeof(double));
        if (!stream.good())
        {
            msg << "Unable to read standard deviation of insert size.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        stream.read(reinterpret_cast<char *>(&histograms[i].readLength), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read read length.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        stream.read(reinterpret_cast<char *>(&histograms[i].offset), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read offset of insert size histogram.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        unsigned histEnd = 0;
        stream.read(reinterpret_cast<char *>(&histEnd), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read insert size histogram size.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        resize(histograms[i].values, histEnd - histograms[i].offset);
        for (unsigned h = 0; h < histEnd - histograms[i].offset; ++h)
        {
            stream.read(reinterpret_cast<char *>(&histograms[i].values[h]), sizeof(double));
            if (!stream.good())
            {
                msg << "Unable to read value from histogram.";
                SEQAN_THROW(ParseError(toCString(msg.str())));
            }
        }
    }

    // Read the number of chromosomes.
    uint32_t numContigs = 0;
    stream.read(reinterpret_cast<char *>(&numContigs), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read number of contigs.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
    resize(contigNames, numContigs);
    resize(contigLengths, numContigs);

    // Read the chromosome names and lengths.
    for (unsigned i = 0; i < numContigs; ++i)
    {
        // Read length of chromosome name.
        uint32_t nameLen = 0;
        stream.read(reinterpret_cast<char *>(&nameLen), sizeof(uint32_t));
        if (!stream.good())
        {
            msg << "Unable to read length of contig name.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }

        // Read the chromosome name.
        resize(contigNames[i], nameLen-1);
        stream.read(reinterpret_cast<char *>(&contigNames[i][0]), nameLen-1);
        stream.read(&buffer[0], 1);
        if (!stream.good())
        {
            msg << "Unable to read contig name.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        if (buffer[0] != '\0') // Expect chromosome names to be null-terminated.
        {
            msg << "Expecting contig name to be null-terminated.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }

        // Read the chromosome length.
        stream.read(reinterpret_cast<char *>(&contigLengths[i]), sizeof(int32_t));
        if (!stream.good())
        {
            msg << "Unable to read contig length.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
    }
}

// =======================================================================================
// Function jumpToRegion()
// =======================================================================================
template<typename TStream>
inline void jumpToRegion(TStream & stream,
                         const String<CharString> & contigNames,
                         const String<int32_t> & contigLengths,
                         unsigned indexRegionSize,
                         const GenomicRegion & region)
{
    uint64_t pos = 0;

    // Add index positions for other chromosomes before the region.
    for (unsigned i = 0; i < length(contigNames); ++i)
    {
        SEQAN_ASSERT_LT(i, length(contigNames));    // Tested before by checkRois().
        if (region.seqName == contigNames[i])
        {
            break;
        }
        else
        {
            pos += contigLengths[i] / indexRegionSize + 1;
        }
    }
    // Add index positions within the region's chromosome.
    pos += region.beginPos / indexRegionSize;

    // Go to the position of the region's file offset in index.
    uint64_t offset = 0;
    stream.seekg(7 + 2 * sizeof(uint32_t) + pos * sizeof(uint64_t));

    // Read the offset and move the stream.
    stream.read(reinterpret_cast<char *>(&offset), sizeof(uint64_t));
    stream.seekg(offset);
}
// ---------------------------------------------------------------------------------------
// Function calculateMedian()
// ---------------------------------------------------------------------------------------
//Calculate the median insert Size from histogram between first and last - 1.
// If the number of counts in hist is even, the lower median is chosen.
inline unsigned calculateMedian(const Iterator<const String<double>, Rooted>::Type & first,
                                const Iterator<const String<double>, Rooted>::Type & last)
{
    double sum = 0;
    Iterator<const String<double>, Rooted>::Type it = first;
    while (it != last)                                         // Sum up all counts in hist
    {
        sum += *it;
        ++it;
    }
    double middle = sum / 2.0;
    it = first;
    double halfSum = 0;
    while (halfSum < middle)                                   // Find index of value separating upper and lower half
    {
        halfSum += *it;
        ++it;
    }
    return position(it) - 1;       // -1 necessary since we stepped over the median with the last increment.
}
// Wrapper for calculating the median on the whole string.
inline unsigned calculateMedian(const String<double> & values)
{
    return calculateMedian(begin(values, Rooted()), end(values, Rooted()));
}
// Wrapper for directly assigning the calculated value to hist.median.
inline void calculateMedian(Histogram & hist,
                            const Iterator<const String<double>, Rooted>::Type & first,
                            const Iterator<const String<double>, Rooted>::Type & last)
{
    hist.median = calculateMedian(first, last);
}
// Wrapper for calculation on whole histogram (except first and last index, which are for irregular values).
inline void calculateMedian(Histogram & hist)
{
    calculateMedian(hist, begin(hist.values, Rooted()) + 1, end(hist.values, Rooted()) - 1);
}
// ---------------------------------------------------------------------------------------
// Function calculateQuantiles()
// ---------------------------------------------------------------------------------------
// Calculates the insert sizes at the 1% and 99% quantile, adjusted for the median.
void calculateQuantiles(Histogram & hist)       // TODO: Will be more efficient if integrated in calculateMedian().
{
    Iterator<const String<double>, Rooted>::Type it = begin(hist.values, Rooted()) + 1;
    Iterator<String<double>, Rooted >::Type itEnd = end(hist.values, Rooted()) - 1;
    unsigned totalCount = 0;
    while (it != itEnd)
    {
        totalCount += *it;
        ++it;
    }
    goBegin(it);
    ++it;
    double lower = totalCount * 0.01;
    double upper = totalCount * 0.99;
    double sum = 0;
    while (it != itEnd)
    {
        sum += *it;
        if (sum >= lower)
        {
            hist.lowerQuantileDist = std::abs(static_cast<int>(position(it)) +
                                     hist.offset - static_cast<int>(hist.median));
            break;
        }
        ++it;
    }
    while (it != itEnd)
    {
        sum += *it;
        if (sum >= upper)
        {
            hist.upperQuantileDist = std::abs(static_cast<int>(position(it)) +
                                     hist.offset - static_cast<int>(hist.median));
            break;
        }
        ++it;
    }
}
// ---------------------------------------------------------------------------------------
// Function _calculateMean()
// ---------------------------------------------------------------------------------------
//Calculate the mean insert size from histogram without considering the offset.
inline void calculateMean(Histogram & hist, THistIter first, const THistIter & last)
{
    double sum = 0;
    double counts = 0;
    while (first != last)
    {
        sum += (*first) * position(first);
        counts += *first;
        ++first;
    }
    if (counts == 0)
    {
        std::ostringstream msg;
        msg << "[PopDel] Error: Attempt to divide by 0 during calculation of histogram mean. Make sure that every read group has sufficient reads in the sampling intervals!";
        SEQAN_THROW(std::range_error(toCString(msg.str())));
    }
    hist.mean = sum/counts;
}
//Wrapper for calculation on whole histogram, except first and last index, which are for irregular values.
inline void calculateMean(Histogram & hist)
{
    calculateMean(hist, begin(hist.values) + 1, end(hist.values) - 1);
}
// ---------------------------------------------------------------------------------------
// Function calculateStddev()
// ---------------------------------------------------------------------------------------
//Calculate the corrected sample standard deviation of the insert size distribution from first to last -1.
// hist.mean has to be calculated before.
inline void calculateStddev(Histogram & histogram, THistIter first, const THistIter & last)
{
    double sum = 0;
    double counts = 0;
    while (first != last)
    {
        sum += (*first) * pow(position(first) - histogram.mean, 2.0);
        counts += *first;
        ++first;
    }
    histogram.stddev = sqrt(sum/(counts - 1));
}
// Wrapper calling the function for the whole histogram (except first and last index, which are for irregular values.)
inline void calculateStddev(Histogram & hist)
{
    calculateStddev(hist, begin(hist.values) + 1, end(hist.values) - 1);
}
inline void calculateHistMetrics(Histogram & hist)
{
    THistIter first = begin(hist.values);
    THistIter last = end(hist.values);
    unsigned minValue = 1;
    unsigned m = length(hist.values) - 1;
    unsigned maxValue = m;
    while (minValue != position(first) || maxValue != position(last))
    {
        first = begin(hist.values) + minValue;
        last = begin(hist.values) + maxValue;
        calculateMean(hist, first, last);
        calculateStddev(hist, first, last);
        calculateMedian(hist, first, last);
        minValue = getHistLeftBorder(hist);
        maxValue = getHistRightBorder(hist);
    }
}
// -----------------------------------------------------------------------------
// Function normalizeValue()
// -----------------------------------------------------------------------------
// Normalize a single histogram value given the insert size window size and read length.
inline void normalizeValue(double & value,
                           const unsigned & insertSize,
                           const unsigned & windowSize,
                           const unsigned & readLength)
{
    int innerDistance = insertSize - 2 * readLength;
    if (innerDistance < 1)
        innerDistance = 1;
    value *=  static_cast<double>(windowSize + innerDistance - 1) / windowSize;
}
// -----------------------------------------------------------------------------
// Function normalizeValues()
// -----------------------------------------------------------------------------
// Normalize all values of a string except the first an last one (=irregular values)
// by calling normalizeValue for each insertSize.
inline void normalizeValues(String<double> & values, const Histogram & hist)
{
    SEQAN_ASSERT_GEQ(hist.offset, 0);
    unsigned insertSize = hist.offset;
    THistIter it = begin(values, Rooted()) + 1;
    THistIter itEnd = end(values, Rooted()) - 1;
    while (it != itEnd)
    {
        normalizeValue(*it, insertSize, hist.windowSize, hist.readLength);
        ++insertSize;
        ++it;
    }
}
// -----------------------------------------------------------------------------
// Function normalizeHistogram()
// -----------------------------------------------------------------------------
//Normalize the histogram, except the irregular values.
inline void normalizeHistogram(Histogram & hist)
{
    normalizeValues(hist.values, hist);             // normalization for single window
}
// -----------------------------------------------------------------------------
// Function smoothHistogram()
// -----------------------------------------------------------------------------
// Smooth the counts in hist by replacing each value with the weighted average of the bases (distance dependend)
// in a window +- 20 BP from its position.
inline void smoothHistogram(Histogram & hist)                      // TODO: Add test
{
    int len = length(hist.values);
    String<double> values;
    resize(values, len);
    for (int i = 0; i < len; ++i)
    {
        double weightedSum = 0;
        double sum = 0;
        for (int j = - 20; j <= 20; ++j)                            // Window of 40 BP (-20, +20)
        {
            double k = std::exp(- j * j / 40.0);                    // Lower weight for distant bases.
            if (i + j >= 0 && i + j < len)
                weightedSum += k * hist.values[i + j];
            sum += k;
        }
        values[i] = weightedSum / sum;
    }
    hist.values = values;
}
// -----------------------------------------------------------------------------
// Function setMinimumProbability()
// -----------------------------------------------------------------------------
// Set minProb to a fraction of the biggest count in values.
// Replace any number in values with minProb if the original value is smaller.
inline void setMinimumProbability(String<double> & values, const double & minProb)
{
    typedef Iterator<String<double>, Standard>::Type TdoubleStringIter;
    for (TdoubleStringIter it = begin(values, Standard()); it < end(values, Standard()); ++it)
        if (*it < minProb)
            *it = minProb;
}
// Set hist.min_prob to a fraction of the biggest count in hist.values.
// Replace any number in hist.values with hist.min_prob if the original value is smaller.
inline void setMinimumProbability(Histogram & hist, const unsigned & pseudoCountFraction)
{
    typedef Iterator<String<double> >::Type TdoubleStringIter;
    double maxValue = 0;
    for (TdoubleStringIter it = begin(hist.values); it < end(hist.values); ++it)  // Get maximum count.
    {
        if (*it > maxValue)
            maxValue = *it;
    }
    hist.min_prob = maxValue / pseudoCountFraction;          //Set min_prob
    for (TdoubleStringIter it = begin(hist.values); it < end(hist.values); ++it) // Replace any value below min_prob.
    {
        if (*it < hist.min_prob)
            *it = hist.min_prob;
    }
}

// -----------------------------------------------------------------------------
// Function printHistogram()
// -----------------------------------------------------------------------------
// Print the Offset and counts of different insert Sizes of histogram to standard output.
inline void printHistogram(const Histogram & hist)
{
    std::cout << "OFFSET=" << hist.offset;
    for (unsigned i = 0; i < length(hist.values); ++i)
        std::cout << "\t" << hist.values[i];
    std::cout << std::endl;
}

// -----------------------------------------------------------------------------
// Function printHistograms()
// -----------------------------------------------------------------------------
// Print insert sizes histograms in human readable format.
template <typename TStream>
void printHistograms(TStream & stream, String<CharString> & readGroups, String<Histogram> & histograms)
{
    // Print the information on insert size distribution.
    stream << "#RG";
    for (unsigned i = 0; i < length(readGroups); ++i)
        stream << "\t" << readGroups[i];
    stream << std::endl;

    stream << "#Median";
    for (unsigned i = 0; i < length(histograms); ++i)
        stream << "\t" << histograms[i].median;
    stream << std::endl;

    stream << "#Std_dev";
    for (unsigned i = 0; i < length(histograms); ++i)
        stream << "\t" << histograms[i].stddev;
    stream << std::endl;

    stream << "#Read_length";
    for (unsigned i = 0; i < length(histograms); ++i)
        stream << "\t" << histograms[i].readLength;
    stream << std::endl;

    // Print the insert size histograms.
    stream << "INSERT_SIZE";
    for (unsigned i = 0; i < length(readGroups); ++i)
        stream << "\t" << readGroups[i];
    stream << std::endl;

    int first = maxValue<int>();
    int last = 0;
    for (unsigned i = 0; i < length(histograms); ++i)
    {
        if (histograms[i].offset < first)
            first = histograms[i].offset;
        if (last < histograms[i].offset + (int)length(histograms[i].values))
            last = histograms[i].offset + length(histograms[i].values);
    }

    for (int i = first; i < last; ++i)
    {
        stream << i;
        for (unsigned rg = 0; rg < length(histograms); ++rg)
        {
            if (i >= histograms[rg].offset && i < histograms[rg].offset + (int)length(histograms[rg].values))
                stream << "\t" << histograms[rg].values[i - histograms[rg].offset];
            else
                stream << "\t" << "NA";
        }
        stream << std::endl;
    }
}

// -----------------------------------------------------------------------------
// Function densityScale()
// -----------------------------------------------------------------------------
// Turn a histogram of counts into a probability density distribution by dividing all values by the sum of all values.
inline void densityScale(String<double> & counts)
{
    unsigned total = 0;
    for (Iterator<String<double>, Standard>::Type it = begin(counts, Standard()); it != end(counts); ++it)
        total += *it;
    for (Iterator<String<double>, Standard>::Type it = begin(counts, Standard()); it != end(counts); ++it)
        *it /= total;
}
// -----------------------------------------------------------------------------
// Function checkUniqueRG()
// -----------------------------------------------------------------------------
// Check if the readGroup has not occured before (= is not in readGroups). Also update the rg counter.
// Return true, if the RG is unique (until now), false otherwise.
inline bool checkUniqueRG(std::map<CharString, unsigned> & readGroups,
                          unsigned & rg,
                          const CharString & readGroup,
                          const CharString & filename)
{
    if (readGroups.count(readGroup) != 0)
    {
        std::ostringstream msg;
        msg << "WARNING: Duplicate insert size histogram of read group \'" 
            << readGroup << "\'. Skipping it in file \'" << filename << "\'.";
        printStatus(msg);
        return false;
    }
    else
    {
        readGroups[readGroup] = rg;
        ++rg;
        return true;
    }
}
// -----------------------------------------------------------------------------
// Function tryOpenHistogram()
// -----------------------------------------------------------------------------
// Trie to open the histogram file and throw an error on failures.
// Return the opened ifstream.
inline std::ifstream tryOpenHistogram(const CharString & filename)
{
    std::ifstream infile(toCString(filename));
    if (!infile.is_open())
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not open histogram file \'" << filename << "\' for reading.";
        SEQAN_THROW(IOError(toCString(msg.str())));
    }
    return infile;
}
// -----------------------------------------------------------------------------
// Function readHistogramLine()
// -----------------------------------------------------------------------------
// First partially clear the histogram hist and then try to read one line of the histogram file.
// The loaded values are written to hist.
// Return 1 on success, -1 if getline() failed, 0 if the RG is a duplicate.
inline int readHistogramLine(Histogram & hist,
                             std::map<CharString, unsigned> & readGroups,
                             unsigned & rg,
                             std::ifstream & infile,
                             const CharString & filename)
{
    hist.clearStrings();
    std::string str;
    if (getline(infile, str))
    {
        std::istringstream iss(str);
        std::string readGroup;
        unsigned last;
        iss >> readGroup;
        iss >> hist.offset;
        iss >> last;
        iss >> hist.median;
        iss >> hist.stddev;
        iss >> hist.readLength;
        if (!checkUniqueRG(readGroups, rg, readGroup, filename))
            return 0;
        double value;
        while(iss >> value)                         // Read the histogram value fields.
            appendValue(hist.values, value);
        return 1;
    }
    else
    {
        return -1;
    }
}
// -----------------------------------------------------------------------------
// Function processHistogram()
// -----------------------------------------------------------------------------
// Prepare the histogram and calculates its metrics.
inline void processHistogram(Histogram & hist,
                             const unsigned & windowSize,
                             const bool & smoothing,
                             const unsigned & pseudoCountFraction)
{
    hist.windowSize = windowSize;
    if (smoothing)
       smoothHistogram(hist);
    calculateQuantiles(hist);
    densityScale(hist.values);
    normalizeHistogram(hist);
    setMinimumProbability(hist, pseudoCountFraction);
}
// -----------------------------------------------------------------------------
// Function loadInsertSizeHistograms()
// -----------------------------------------------------------------------------
// Load histogram of insert sizes from file and call _normalizeHistogram, _setMinimumProbability, _logHistogram.
// Add the loaded and processed histograms to the set of all histograms (one hist per Read Group)
inline void loadInsertSizeHistograms(String<Histogram> & histograms,           // TODO: Add test. Rename to loadHeaders?
                                     std::map<CharString, unsigned> & readGroups,
                                     String<unsigned> & rgs,
                                     const CharString & filename,
                                     String<String<CharString> > & contigNames,
                                     String<String<int32_t> > & contigLengths,
                                     String<unsigned> & indexRegionSizes,
                                     bool smoothing,
                                     const unsigned & pseudoCountFraction)
{
    std::ifstream infile = tryOpenHistogram(filename);
    String<CharString> sampleReadGroups;
    String<Histogram> sampleHistograms;
    String<CharString> sampleContigNames;
    String<int32_t> sampleContigLengths;  // TODO: Make this more efficient by directly writing to the final sets.
    unsigned indexRegionSize = 0;
    readProfileHeader(infile,
                      filename,
                      sampleReadGroups,
                      sampleHistograms,
                      sampleContigNames,
                      sampleContigLengths,
                      indexRegionSize);
    // Process all histograms of the sample and add them.
    unsigned rg = length(readGroups);
    for (unsigned i = 0; i < length(sampleReadGroups); ++i)
    {
        if (checkUniqueRG(readGroups, rg, sampleReadGroups[i], filename))
        {
            appendValue(rgs, rg - 1);
            processHistogram(sampleHistograms[i], 256, smoothing, pseudoCountFraction);
            appendValue(histograms, sampleHistograms[i]);
        }
    }
    appendValue(indexRegionSizes, indexRegionSize);
    appendValue(contigNames, sampleContigNames);
    appendValue(contigLengths, sampleContigLengths);
}
//Overload for applying loadInsertSizeHistograms on multiple files.
inline void loadInsertSizeHistograms(String<Histogram> & histograms,
                                     std::map<CharString, unsigned> & readGroups,
                                     TRGs & rgs,
                                     const String<CharString> & filenames,
                                     String<String<CharString> > & contigNames,
                                     String<String<int32_t> > & contigLengths,
                                     String<unsigned> & indexRegionSizes,
                                     bool smoothing,
                                     const unsigned & pseudoCountFraction)
{
    resize(rgs, length(filenames), Exact());
    unsigned sampleNum = 0;
    for (unsigned i = 0; i < length(filenames); ++i)
    {
        loadInsertSizeHistograms(histograms,
                                 readGroups,
                                 rgs[sampleNum],
                                 filenames[i],
                                 contigNames,
                                 contigLengths,
                                 indexRegionSizes,
                                 smoothing,
                                 pseudoCountFraction);
        std::ostringstream msg;
        msg << "Loaded histogram from file \'" << filenames[i] << "\'.";
        printStatus(msg);
        ++sampleNum;
    }
}
// -----------------------------------------------------------------------------
// Function I()
// -----------------------------------------------------------------------------
// Take the insert size deviation and
// return the value stored in hist.values[i] if i lies between min and max insert size in the histogram.
// Return hist.min_prob otherwise.
inline double I(const Histogram & hist, int deviation)
{
    int i = deviation + static_cast<int>(hist.median) - hist.offset;
    if (i <= 0 || i + 1 >= static_cast<int>(length(hist.values)))       // Exlclude first and last index.
        return (hist.min_prob);
    return hist.values[i];
}
// -----------------------------------------------------------------------------
// Function insertSizeToHistValue()
// -----------------------------------------------------------------------------
// Take the insert size (not the deviation!) and
// return the value stored in hist.values[i] if i lies between min and max insert size in the histogram.
// Return hist.min_prob otherwise.
inline double insertSizeToHistValue(const Histogram & hist, int InsertSize)
{
    int i = InsertSize - hist.offset;
    if (i <= 0 || i + 1 >= static_cast<int>(length(hist.values)))       // Exlclude first and last index.
        return (hist.min_prob);
    return hist.values[i];
}
#endif /* INSERT_HISTOGRAM_POPDEL_H_ */
