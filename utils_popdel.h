#ifndef UTILS_POPDEL_H_
#define UTILS_POPDEL_H_

#include<ctime>
#include <unordered_set>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

using namespace seqan;

// =======================================================================================
// Typedefs
// =======================================================================================
typedef __uint32 TEntryIdType;                            // Type for id in ChromosomeProfileEndEntry;
typedef String<int> TProfile;                             // Profile values of a single RG in a single window
typedef String<TProfile> TReadGroupProfiles;              // Profiles per ReadGroups in a single Window
typedef String<TReadGroupProfiles> TWindowProfiles;       // 3D Array: Window, Readgroup, ProfileValues
typedef Iterator<TWindowProfiles>::Type TWindowIter;
typedef String<TEntryIdType> TReadGroupIndices;
typedef String<TReadGroupIndices> TRGs;
// =======================================================================================
// Classes
// =======================================================================================
// =======================================================================================
// Class Dad
// =======================================================================================
// Class for counting how the insertSizes of a read group overlap the different histograms.
struct Dad
{
    unsigned ref;       // only in ref hist (or left of it).
    unsigned both;      // in both hists.
    unsigned between;   // between both hists, but in neither of them.
    unsigned alt;       // only in alternative hist.
    unsigned right;     // right of alternative hist.

    Dad() : ref(0), both(0), between(0), alt(0), right(0){}

    Dad(unsigned re, unsigned bo, unsigned be, unsigned al, unsigned ri) :
    ref(re),
    both(bo),
    between(be), 
    alt(al),
    right(ri){}

    inline void reset()
    {
        ref = both = between = alt = right = 0u;
    }
};
// Object storing all information associated with a deletion call
// =======================================================================================
// Class Call
// =======================================================================================
// Object storing all information associated with a deletion call
struct Call
{
    __uint32 initialLength;         // Estimate for deletion length at initialization.
    __uint32 iterations;            // Number of iterations performed.
    __uint32 deletionLength;        // Final estimate for deletion length.
    double   likelihoodRatio;       // Likelihood ratio calculated for this deletion.
    String<Triple<unsigned> > lads; // Likelihood based counts of pairs supporting: i1: Ref, i2: Ambiguous, i3: Del. One per Sample.
    String<Dad> dads; // distribution based counts of pairs supporting:
                                    //  0: Ref, 1: both, 2: between but none, 3: Del., 4: bigger than del. One per Sample.
    String<Pair<unsigned> > firstLast;// Pos of the lowest first and the highest last win of all active reads.
    double   frequency;             // Allele-frequency of the deletion across all samples.
    unsigned position;              // Assumed start-position (window) of the variant.
    unsigned endPosition;           // Assumed end-position (window) of the variant.
    String<Triple<unsigned> > gtLikelihoods; // PHRED-scaled GT likelihhods of each sample. Order: HomRef, Het, HomDel.
    unsigned char filter;           // 8 Bits indicating the failed filters. 1 at the position indicates failed filter.
                                    // (right to left) 1.: LR ratio test failed, 2.: High coverage, 3.:Sample Number
    Call() :
    initialLength(0),
    iterations(0),
    deletionLength(0),
    likelihoodRatio(0.0),
    frequency(0.0),
    position(0),
    endPosition(0),
    filter(0){}

    Call(__uint32 initLen,
         __uint32 itNum,
         __uint32 delLen,
         double llr,
         double f,
         unsigned pos,
         unsigned endPos = 0) :
    initialLength(initLen),
    iterations(itNum),
    deletionLength(delLen),
    likelihoodRatio(llr),
    frequency(f),
    position(pos),
    endPosition(endPos){}

    void reset(void)
    {
        initialLength = 0;
        iterations = 0;
        deletionLength = 0;
        likelihoodRatio = 0;
        clear(lads);
        clear(dads);
        frequency = 0;
        position = 0;
        endPosition = 0;
        filter = 0;
    }
};
// =======================================================================================
// Functions
// =======================================================================================
// =======================================================================================
// Function setLRFilter()
// =======================================================================================
// Set the bit for the likelihood ratio filter.
inline void setLRFilter(Call & call)
{
    call.filter |= 1;   // 0001
}
// =======================================================================================
// Function checkLRPass()
// =======================================================================================
// Return true if the likelihood ratio filter has been passed, false otherwise.
inline bool checkLRPass(const Call & call)
{
   return !(call.filter & 1);     // 0001
}
// =======================================================================================
// Function setCoverageFilter()
// =======================================================================================
// Set the bit for the high coverage filter.
inline void setCoverageFilter(Call & call)
{
    call.filter |= 2;   // 0010
}
// =======================================================================================
// Function checkCoveragePass()
// =======================================================================================
// Return true if the high coverage filter has been passed, false otherwise.
inline bool checkCoveragePass(const Call & call)
{
    return !(call.filter & 2);   // 0010
}
// =======================================================================================
// Function setSampleFilter()
// =======================================================================================
// Set the bit for the sample number filter.
inline void setSampleFilter(Call & call)
{
    call.filter |= 4;   // 0100
}
// =======================================================================================
// Function checkSamplePass()
// =======================================================================================
// Return true if the sample number filter has been passed, false otherwise.
inline bool checkSamplePass(const Call & call)
{
    return !(call.filter & 4);   // 0100
}
// =======================================================================================
// Function setGT0Filter()
// =======================================================================================
// Set the bit for the sample number filter.
inline void setGT0Filter(Call & call)
{
    call.filter |= 8;   // 1000
}
// =======================================================================================
// Function checkGT0Pass()
// =======================================================================================
// Return true if the sample number filter has been passed, false otherwise.
inline bool checkGT0Pass(const Call & call)
{
    return !(call.filter & 8);   // 1000
}
// =======================================================================================
// Function checkAllPass()
// =======================================================================================
// Return true if the high coverage filter has been passed, false otherwise.
inline bool checkAllPass(const Call & call)
{
    return checkCoveragePass(call) & checkLRPass(call) & checkSamplePass(call) && checkGT0Pass(call);
}
// Reset all filters.
inline void resetFilters(Call & call)
{
    call.filter = 0;
}
// =======================================================================================
// Function markInvalidCall()
// =======================================================================================
// Marks a call as invalid by setting its filter to maxValue<unsigned>().
inline void markInvalidCall(Call & call)
{
    call.filter = maxValue<unsigned char>();
}
// =======================================================================================
// Function checkDelSizeSimilar()
// =======================================================================================
// Return true if the two delSizes a and b are within f percent of each other or within +-2 sddev.
inline bool checkDelSizeSimilar(const unsigned & a, const unsigned & b, const double & stddev, const double f = 0.8)
{
    unsigned l;
    unsigned r;
    if (a < b)
    {
        l = a;
        r = b;
    }
    else
    {
        l = b;
        r = a;
    }
    if (l + 2 * stddev >= r)
        return true;
    return l >= f * r;
}
// =======================================================================================
// Function checkEnoughOverlap()
// =======================================================================================
// Return true if the calls a and b overlap for at least for f-percent of the smaller variant's number of windows.
inline bool checkEnoughOverlap(const Call & a, const Call & b, const double f = 0.5)
{   //TODO: Consider making this more generous for smaller deletions
    unsigned minLen = std::min(a.deletionLength, b.deletionLength);
    unsigned left = std::max(a.position, b.position);
    unsigned right = std::min(a.position + a.deletionLength, b.position + b.deletionLength);
    int overlap = (right - left);
    return overlap >= f * minLen;
}
// =======================================================================================
// Function areSimilar()
// =======================================================================================
// Return true, if the two given SupportStretch are similar enough to be considered equal.
inline bool similar(const Call & a, const Call & b, const double & stddev)
{
    return checkDelSizeSimilar(a.deletionLength, b.deletionLength, stddev) && checkEnoughOverlap(a, b);
}
// =======================================================================================
// Function lowerCall()
// =======================================================================================
// Return true if the first Call comes not after the second, false otherwise.
// The First the positions are compared and for equal positions the shorter variant comes before the longer.
// If both are equal, the element with the higher LR is prefered.
inline bool lowerCall(const Call & l, const Call & r)
{
    if (l.position < r.position)
        return true;
    else if (l.position > r.position)
        return false;
    else if (l.deletionLength < r.deletionLength)
        return true;
    else if (l.deletionLength > r.deletionLength)
        return false;
    else if (l.likelihoodRatio > r.likelihoodRatio)
        return true;
    else
        return false;
}
// =======================================================================================
// Function unifyCalls()
// =======================================================================================
// Unifies duplicate calls in c.
// Return false if the string of calls is empty, true otherwise.
inline bool unifyCalls(String<Call> & calls, const double & stddev)
{
    if (length(calls) == 0)
        return false;
    std::sort(begin(calls), end(calls), lowerCall);
    unsigned callCount = 1;

    Iterator<String<Call> >::Type currentIt = begin(calls, Standard());
    String<Triple<unsigned> > genotypes;
    resize(genotypes, length(currentIt->gtLikelihoods), Triple<unsigned>(0, 0, 0));
    unsigned winCount = 1;
    for (Iterator<String<Call>, Standard >::Type it = begin(calls) + 1; it != end(calls); ++it)
    {
        if (similar(*currentIt, *it, stddev))
        {
            // Count genotypes for all samples across all windows and take their means.
            ++winCount;
            for (unsigned i = 0; i < length(genotypes); ++i)
            {
                genotypes[i].i1 += it->gtLikelihoods[i].i1;
                genotypes[i].i2 += it->gtLikelihoods[i].i2;
                genotypes[i].i3 += it->gtLikelihoods[i].i3;
            }
            markInvalidCall(*it);
        }
        else if (winCount != 1)
        {
            for (unsigned i = 0; i < length(genotypes); ++i)
            {
                double minGt = std::min(std::min(genotypes[i].i1, genotypes[i].i2), genotypes[i].i3);
                currentIt->gtLikelihoods[i].i1 = std::round(static_cast<double> (genotypes[i].i1 - minGt) / winCount);
                currentIt->gtLikelihoods[i].i2 = std::round(static_cast<double> (genotypes[i].i2 - minGt) / winCount);
                currentIt->gtLikelihoods[i].i3 = std::round(static_cast<double> (genotypes[i].i3 - minGt) / winCount);
                genotypes[i].i1 = 0;
                genotypes[i].i2 = 0;
                genotypes[i].i3 = 0;
            }
            currentIt = it;
            ++callCount;
            winCount = 1;
        }
    }
    String<Call> tmp;
    reserve(tmp, callCount, Exact());
    for (Iterator<String<Call> >::Type it = begin(calls); it != end(calls); ++it)
    {
        if (it->filter != maxValue<unsigned char>())
            appendValue(tmp, *it);
    }
    move(calls, tmp);
    return true;
}
// =======================================================================================
// Function sum()
// =======================================================================================
// Return the sum of the two values in a pair of doubles.
inline double sum(const Pair<double> & p)
{
    return (p.i1 + p.i2);
}
// Return the sum of all three values in a triple of unsigned integers.
inline unsigned sum(const Triple<unsigned> & t)
{
    return (t.i1 + t.i2 + t.i3);
}
// =======================================================================================
// Function allTrue()
// =======================================================================================
// Return true, if all bools in the string are true.
inline bool allTrue(const String<bool> & s)
{
    for (Iterator<const String<bool> >::Type it = begin(s, Standard()); it != end(s); ++it)
    {
        if (!*it)
            return false;
    }
    return true;
}
// =======================================================================================
// Function checkSampleNumber()
// =======================================================================================
// Checks if the number of samples with data at the current position is high enough.
// Return true if so, false otherwise.
inline bool checkSampleNumber(const Call & call, const double & threshold)
{
    unsigned total = length(call.gtLikelihoods);
    unsigned dataSample = total;
    for (unsigned  i = 0; i < total; ++i)
        if (sum(call.gtLikelihoods[i]) == 0)
            --dataSample;
    if (static_cast<double>(dataSample) / total < threshold)
        return false;
    else
        return true;
}
inline int32_t max(const String<int32_t> & s)
{
    int32_t currentMax = 0;
    for (Iterator<const String<int32_t>, Standard>::Type it = begin(s, Standard()); it != end(s); ++it)
    {
        if (*it > currentMax)
            currentMax = *it;
    }
    return currentMax;
}
// =======================================================================================
// Function printStatus()
// =======================================================================================
// Print the status message plus date and time.
void printStatus(const char * message)
{
    time_t now = time(0);                                                       //Get the current date and time.
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "[PopDel %Y-%m-%d %X] ", &tstruct);
    std::cout << buf << message << std::endl;                                   //Print time and message.
}
// Wrapper for calling printStatuts with ostringstream instead of CString.
void printStatus(std::ostringstream & message)
{
    std::string msg = message.str();
    printStatus(toCString(msg));
}
// =======================================================================================
// Function loadBai()
// =======================================================================================
// Load the BAI-file belonging to the given BAM-file. 
// Return false on errors, true otherwise.
template<typename TString>
inline bool loadBai(BamIndex<Bai> & bai, TString filename)
{
    filename += ".bai";
    if (!open(bai, toCString(filename)))
    {
        CharString message = "[PopDel] ERROR: Could not load the BAI file \'";
        message += filename;
        message += "\'.";
        SEQAN_THROW(IOError(toCString(message)));
        return false;
    }
    return true;
}
// =======================================================================================
// Function parseGenomicRegion()
// =======================================================================================
// Wrapper aroung SeqAns parse() function to guarantee conformity with samtools style regions.
// Return false on parsing errors, true otherwise.
template<typename TString>
inline bool parseGenomicRegion(GenomicRegion & reg, const TString & s)
{
    if (length(s) == 0)
    {
        std::ostringstream msg;
        msg << "[PopDel] ERROR: Empty genomic region.";
        printStatus(msg);
        return false;
    }
    try
    {
        parse(reg, s);      //Try regular parsing first.
        if (!isalnum(reg.seqName[length(reg.seqName) - 1]))
        {
            std::ostringstream msg;
            msg << "[PopDel] ERROR: Invalid genomic region \'" << s << "\'.";
            printStatus(msg);
            return false;
        }
        return true;
    }
    catch(...)
    {
        unsigned lastCol = 0;                                                         // Store the last occurenc of ':'.
        unsigned i = 0;
        for (i = 0; i < length(s); ++i)
        {
            if (s[i] == ':')
                lastCol = i;
        }
        if (lastCol == i - 1) // last char in s i ':' -> invalid
        {
            std::ostringstream msg;
            msg << "[PopDel] ERROR: Invalid genomic region \'" << s << "\'.";
            printStatus(msg);
            return false;
        }
        if (!isdigit(s[lastCol + 1]))                      // First character after last ':' not a digit -> all seqName.
        {
            if (!isalnum(s[length(s) - 1]))                 // Is last char of s valid?
            {
                std::ostringstream msg;
                msg << "[PopDel] ERROR: Invalid genomic region \'" << s << "\'.";
                printStatus(msg);
                return false;
            }
            reg.seqName = s;
            return true;
        }
        typename Prefix<const  TString>::Type first = prefix(s, lastCol);     // Prefix up to last :
        typename Suffix<const TString>::Type last = suffix(s, lastCol + 1);  // Suffix from last :
        CharString dummy = "x:";
        append(dummy, last);
        try
        {
            parse(reg, dummy);
        }
        catch(...)
        {
            std::ostringstream msg;
            msg << "[PopDel] ERROR: Invalid genomic region \'" << s << "\'.";
            printStatus(msg);
            return false;
        }
        if (!isalnum(first[length(first) - 1]))
        {
            std::ostringstream msg;
            msg << "[PopDel] ERROR: Invalid genomic region \'" << s << "\'.";
            printStatus(msg);
            return false;
        }
        reg.seqName = first;
        return true;
    }
}

// =======================================================================================
// Function loadFilenames()
// =======================================================================================
// Load multiple filenames from a single file. Removes path to file containing the filenames from set.
void loadFilenames(String<CharString> & files)
{
    SEQAN_ASSERT_EQ(length(files), 1u);
    CharString inputfilename = files[0];
    clear(files);
    // Open infile.
    std::ifstream infile(toCString(inputfilename));
    if (!infile.is_open())
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not open file \'" << inputfilename << "\' for reading.";
        SEQAN_THROW(IOError(toCString(msg.str())));
    }
    // Read file names from infile.
    std::string filename;
    std::unordered_set<std::string> tmpSet;
    while (infile >> filename)
    {   //First add them to a set to avoid duplicates.
        if (tmpSet.insert(filename).second)
        {
            appendValue(files, filename);
        }
        else
        {
            std::ostringstream msg;
            msg << "WARNING: duplicate file " << filename << ". Ignoring additional occurences.";
            printStatus(msg);
        }
    }
    // Print status message.
    std::ostringstream msg;
    msg << "Loaded " << length(files) << " filenames from \'" << inputfilename << "\'.";
    printStatus(msg);
}

// =======================================================================================
// Function _getReadGroup()
// =======================================================================================
// Extract the read-group encoded in tags.
// Return the read-group as a CharString.
inline CharString getReadGroup(const CharString & tags)
{
    BamTagsDict dict(tags);
    unsigned key = 0;
    findTagKey(key, dict, "RG");                                     //TODO faster access if position in dict is known
    CharString rg = "";
    extractTagValue(rg, dict, key);
    return rg;
}
// Extract the read-group encoded in tags
// Return its rank in the header by extracting it from in the map of read groups.
inline unsigned getReadGroup(const CharString & tags, const std::map<CharString, unsigned> & readGroups)
{
    CharString rg = getReadGroup(tags);
    SEQAN_ASSERT_NEQ(readGroups.count(rg), 0u);
    return readGroups.at(rg);
}
// =======================================================================================
// Function _getReadGroups()
// =======================================================================================
// Extract all read group IDs and their rank in the header and write them to map.
// Return the number of read groups in the header.
inline unsigned getReadGroups(std::map<CharString, unsigned> & readGroups, const BamHeader & header)
{
    typedef Iterator<const BamHeader>::Type THeaderIter;
    readGroups.clear();
    unsigned k = 0;
    THeaderIter itEnd = end(header);
    for (THeaderIter it = begin(header); it != itEnd; ++it)
    {
        if ((*it).type == BAM_HEADER_READ_GROUP)
        {
            unsigned idx = 0;
            findTagKey(idx, "ID", *it);                        //Find the read-Group ID-tag...
            CharString rg = "";
            getTagValue(rg, idx, *it);                          //... and get the read-group-ID stored in ID-tag
            readGroups[rg] = k;                                 //Write ID and rank of this read group in header to map
            ++k;
        }
    }
    SEQAN_ASSERT_EQ(readGroups.size(), k);
    return k;
}
// =======================================================================================
// Functions for genomic interval loading and processing
// =======================================================================================
// =======================================================================================
// Function _getWholeGenomeIntervals()
// =======================================================================================
// Append start and end positions + sequence names of all contigs/chromosomes and store them in String of GenomicRegion.
// Return total length of all intervals
inline unsigned getWholeGenomeIntervals(String<GenomicRegion> & intervals, const BamHeader & header)
{
    typedef Iterator<const BamHeader>::Type THeaderIter;
    unsigned totalLength = 0;
    THeaderIter itEnd = end(header);
    for (THeaderIter it = begin(header); it != itEnd; ++it)
    {
        if ((*it).type == BAM_HEADER_REFERENCE)                         //Reference Sequence Dictionary (@SQ)
        {
            GenomicRegion itv;
            itv.seqName = (*it).tags[0].i2;                             //Store chromosome/contig name of region
            itv.beginPos = 0;
            lexicalCastWithException(itv.endPos, (*it).tags[1].i2);     //Get the end position of the interval
            appendValue(intervals, itv);                                //Append new interval to list of intervals
            totalLength += itv.endPos;                                  //Increase total length of sum of intervals
        }
    }
    return totalLength;
}
// =======================================================================================
// Struct lowerGenomicRegion()
// =======================================================================================
// Return true if the first GenomicRegion starts before the second one starts.
// Meant for stable_sort.
inline bool lowerGenomicRegion(const GenomicRegion & r1, const GenomicRegion & r2)
{
    Lexical<> cmp(r1.seqName, r2.seqName);                      //Lexicographically compare the sequence names.
    if (isLess(cmp))                                            //r1.seqName is lexicographically smaller
        return true;
    else if (isGreater(cmp))                                    //r1.seqName is lexicographically bigger
        return false;
    else if (r1.beginPos <= r2.beginPos)
        return true;
    else
        return false;
}
// =======================================================================================
// Function _fillInvalidPositions()
// =======================================================================================
// Replace invalid values of GeomicRegions object with processable values.
// Use 0 for start position. Use end position of sequence for end position of interval.
inline void fillInvalidPositions(GenomicRegion & itv, const BamHeader & header)
{
    typedef Iterator<const BamHeader>::Type THeaderIter;
    if (itv.beginPos == GenomicRegion::INVALID_POS)
        itv.beginPos = 0;
    if (itv.endPos == GenomicRegion::INVALID_POS)
    {
        THeaderIter itEnd = end(header);
        for (THeaderIter it = begin(header); it != itEnd; ++it)
        {
            if ((*it).type == BAM_HEADER_REFERENCE && (*it).tags[0].i2 == itv.seqName)
                lexicalCastWithException(itv.endPos, (*it).tags[1].i2);                     //End position of sequence
        }
    }
}

inline Iterator<const String<GenomicRegion> >::Type findInterval(const String<GenomicRegion> & intervals,
                                                                 const GenomicRegion & roi)
{
    for (Iterator<const String<GenomicRegion> > ::Type it = begin(intervals); it != end(intervals); ++it)
    {
        if (it->seqName == roi.seqName)
            return it;
    }
    return end(intervals);
}
// ---------------------------------------------------------------------------------------
// Function createRegularIntervals()
// ---------------------------------------------------------------------------------------
// Return a string of the default sampling regions for profiling.
void createRegularIntervals(String<GenomicRegion> & intervals, const String<GenomicRegion> & rois)
{
    String<GenomicRegion> tmpIntervals;

    resize(tmpIntervals, 22, Exact());
    parse(tmpIntervals[0], "chr1:35000000-36000000");
    parse(tmpIntervals[1], "chr2:174000000-175000000");
    parse(tmpIntervals[2], "chr3:36500000-37500000");
    parse(tmpIntervals[3], "chr4:88000000-89000000");
    parse(tmpIntervals[4], "chr5:38000000-39000000");
    parse(tmpIntervals[5], "chr6:38000000-39000000");
    parse(tmpIntervals[6], "chr7:38000000-39000000");
    parse(tmpIntervals[7], "chr8:19000000-20000000");
    parse(tmpIntervals[8], "chr9:19000000-20000000");
    parse(tmpIntervals[9], "chr10:19000000-20000000");
    parse(tmpIntervals[10], "chr11:19000000-20000000");
    parse(tmpIntervals[11], "chr12:19000000-20000000");
    parse(tmpIntervals[12], "chr13:25000000-26000000");
    parse(tmpIntervals[13], "chr14:25000000-26000000");
    parse(tmpIntervals[14], "chr15:25000000-26000000");
    parse(tmpIntervals[15], "chr16:25000000-26000000");
    parse(tmpIntervals[16], "chr17:31000000-32000000");
    parse(tmpIntervals[17], "chr18:31000000-32000000");
    parse(tmpIntervals[18], "chr19:31000000-32000000");
    parse(tmpIntervals[19], "chr20:33000000-34000000");
    parse(tmpIntervals[20], "chr21:21000000-22000000");
    parse(tmpIntervals[21], "chr22:25000000-26000000");

    for (Iterator<const String<GenomicRegion> >::Type roiIt = begin(rois); roiIt != end(rois); ++roiIt)
    {
        Iterator<const String<GenomicRegion> >::Type itvIt = findInterval(tmpIntervals, *roiIt);
        if (itvIt != end(tmpIntervals))
        {
            if(findInterval(intervals, *roiIt) == end(intervals))
            {   //Only append the interval if it has not already been added.
                appendValue(intervals, *itvIt);
            }
        }
    }

    if (empty(intervals))
    {
        std::ostringstream msg;
        msg << "[PopDel] No contig name of any ROI matches any of the contig names of the default intervals for"
               " the parameter estimation. Please check the contig names of the ROI's and/or use user-defined"
               " sampling regions (option \'-i\')";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
}
// =======================================================================================
// Function readIntervals()
// =======================================================================================
// Load genomic intervals from a file and store them in a string of GenomicRegion objects.
inline void readIntervals(String<GenomicRegion> & intervals,
                          const CharString & filename,
                          const BamHeader & header,
                          const String<GenomicRegion> & rois)
{
    if (filename != "")
    {
        size_t initialLength = length(intervals);
        std::ifstream infile(toCString(filename));
        if (!infile.is_open())
        {
            std::ostringstream msg;
            msg << "[PopDel] Could not open interval file \'" << filename << "\' for reading.";
            SEQAN_THROW(IOError(toCString(msg.str())));
        }
        std::string word;
        GenomicRegion itv;
        while (infile >> word)
        {
            parseGenomicRegion(itv, word);                                                           //rIDs will not be set
            fillInvalidPositions(itv, header);
            appendValue(intervals, itv);
        }
        std::ostringstream msg;
        msg << "Finished reading " << (length(intervals) - initialLength) << " intervals from file \'" << filename << "\'.";
        printStatus(msg);
    }
    else
    {
        createRegularIntervals(intervals, rois);
    }
}
// =======================================================================================
// Function expandIntervals()
// =======================================================================================
// Expands the intervals by the deletionSize.
inline void expandInterals(String<GenomicRegion> & intervals, const int & maxDeletionLength)
{
    Iterator<String<GenomicRegion>, Rooted >::Type it = begin(intervals, Rooted());
    while (!atEnd(it))
    {
        if (it->beginPos > maxDeletionLength)
            it->beginPos -= maxDeletionLength;
        else
            it->beginPos = 0;
        it->endPos += maxDeletionLength;
        ++it;
    }
}

// =======================================================================================
// Function parseIntervals()
// =======================================================================================
// Parse genomic intervals (if given in command line instead of file) and store them in string of GenomicRegion objects.
inline void parseIntervals(String<GenomicRegion> & intervals,
                           const std::vector<std::string> & intervalStrings,
                           const CharString & bamfile)
{
    typedef std::vector<std::string>::const_iterator TIter;
    BamFileIn infile;
    open(infile, toCString(bamfile));
    BamHeader header;
    readHeader(header, infile);
    TIter it = intervalStrings.begin();
    TIter itEnd = intervalStrings.end();
    GenomicRegion itv;
    while (it != itEnd)
    {
        if(!parseGenomicRegion(itv, *it))                                                    //rIDs will not be set
        {
            std::ostringstream msg;
            msg << "[PopDel] Error while parsing genomic region \'" << *it << "\'. Terminating.";
            SEQAN_THROW(ParseError(toCString(msg.str())));
        }
        fillInvalidPositions(itv, header);
        appendValue(intervals, itv);
        ++it;
    }
}
// =======================================================================================
// Function mergeOverlappingIntervals()
// =======================================================================================
// Use GenomicRegionLess for sorting the intervals in the string of GenomicRegion objects and merge the overlapping ones.
inline void mergeOverlappingIntervals(String<GenomicRegion> & intervals, int maxDeletionLength=0)
{
    typedef Iterator<String<GenomicRegion> >::Type TIter;
    if (length(intervals) < 2)
        return;
    std::stable_sort(begin(intervals), end(intervals), &lowerGenomicRegion);
    TIter itEnd = end(intervals);
    TIter it = begin(intervals);
    String<GenomicRegion> mergedIntervals;
    appendValue(mergedIntervals, *it);
    TIter last = begin(mergedIntervals);
    ++it;
    while (it < itEnd)
    {
        if ((*last).seqName == (*it).seqName && (*last).endPos + maxDeletionLength >= ((*it).beginPos))
        {
            std::ostringstream msg;
            msg << "WARNING: Intervals \'" << (*last).seqName;
            if ((*last).beginPos >= 0)
            {
                if((*last).beginPos != 0)
                {
                    msg << ":" << (*last).beginPos + 1;
                }
                if ((*last).endPos != maxValue<__int32>())
                {
                    if((*last).beginPos == 0)
                    {
                        msg << ":" << (*last).beginPos + 1 << "-" << (*last).endPos;
                    }
                    else
                    {
                        msg << "-" << (*last).endPos;
                    }
                }
            }
            msg << "\' and \'" << (*it).seqName;
            if ((*it).beginPos >= 0)
            {
                msg << ":" << (*it).beginPos + 1;
                if ((*it).endPos != maxValue<__int32>())
                {
                    msg << "-" << (*it).endPos;
                }
            }
            msg << "\' overlap or lie within " << maxDeletionLength <<"-bp vicinity and have been merged.";
            if ((*last).endPos < (*it).endPos)
                (*last).endPos = (*it).endPos;
            printStatus(msg);
        }
        else
        {
            appendValue(mergedIntervals, *it);
            ++last;
        }
        ++it;
    }
    SEQAN_ASSERT_LEQ(length(mergedIntervals), length(intervals));
    intervals = mergedIntervals;
}
// =======================================================================================
// Function _setRIDs()
// =======================================================================================
// Check if the seqName of each interval is also present in the BAM-file. If so, set the rID of the interval.
// If not, set rID to INVALID_REFID.
inline void setRIDs(String<GenomicRegion> & intervals, const BamFileIn & infile)
{                                                   //TODO: Check if regions without valid ID is realy skipped lateron
    typedef Iterator<String<GenomicRegion> >::Type TItvIter;
    TItvIter itv = begin(intervals);
    TItvIter itvEnd = end(intervals);
    while (itv != itvEnd)
    {
        if (!getIdByName((*itv).rID, contigNamesCache(context(infile)), (*itv).seqName))
        {
            (*itv).rID = BamAlignmentRecord::INVALID_REFID;
            CharString itvStr;
            (*itv).toString(itvStr);
            std::cerr << "WARNING: Chromosome/contig name unknown. Skipping interval \'";
            std::cerr << itvStr << "\'." << std::endl;
        }
        ++itv;
    }
}
// =======================================================================================
// Function _round()
// =======================================================================================
// Round d mathematically to the nearest integer.
inline int _round(double d)
{
    return std::floor(d + 0.5);
}
// =======================================================================================
// Function getContigRank()
// =======================================================================================
// Takes a contig name and returns it's relative rank (i.e. distance to allROIsBegin) in the string of all ROIs.
// Checks only the contig names between the two iterators (excluding the last position.)
// If contigName is not present, maxValue<unsigned>() is returned.
inline unsigned getContigRank(const CharString & contigName,
                              Iterator<const String<GenomicRegion> >::Type allROIsBegin,
                              const Iterator<const String<GenomicRegion> >::Type & allROIsEnd)
{
    unsigned i = 0;
    while (allROIsBegin != allROIsEnd)
    {
        if (contigName == allROIsBegin->seqName)
            return i;
        else
            ++allROIsBegin;
        ++i;
    }
    return maxValue<unsigned>();
}
// =======================================================================================
// Function lowerCoord()
// =======================================================================================
// Take 2 coordinates (chrom, pos) and return true if the first one comes before the second one (or if they are equal).
// Use nextROI and allRois to determine the rank of the contig.
inline bool lowerCoord(const Pair<CharString, unsigned> & first,
                       const Pair<CharString, unsigned> & second,
                       const Iterator<String<GenomicRegion>, Standard>::Type & nextROI,
                       const String<GenomicRegion> & allRois)
{
    const Iterator<const String<GenomicRegion>, Standard>::Type ROIend = end(allRois, Standard());
    unsigned rankFirst = getContigRank(first.i1, nextROI, ROIend);
    unsigned rankSecond = getContigRank(second.i1, nextROI, ROIend);
    Lexical<> cmp(first, second);                      //Lexicographically compare the sequence names.
    if (rankFirst < rankSecond)                                   //first is lexicographically smaller.
        return true;
    else if (rankFirst > rankSecond)                           //first is lexicographically greater.
        return false;
    else if (first.i2 <= second.i2)
        return true;
    else
        return false;
}
// =======================================================================================
// Function isValidCall()
// =======================================================================================
// Check if a call in the buffer is valid or not. This is necessary since the buffer gets not completely cleared.
inline bool isValidCall(const Call & call)
{
    return call.iterations != 0;        // only true if the call has been changed since the last pass.
}
// =======================================================================================
// Function getSampleName()
// =======================================================================================
// Extract the last part of the path to get the file/sample name without file format ending.
inline Infix<const String<char> >::Type getSampleName(const String<char> & path)
{
    String<char> filename;
    unsigned lastDelimPos = 0;
    unsigned lastDotPos = length(path);
    Iterator<const String<char>, Rooted>::Type it = begin(path, Rooted());
    while (!atEnd(it))
    {
        if (*it == '/' || *it == '\\')
            lastDelimPos = position(it);               // Save last occurence of delimiter
        else if (*it == '.')
            lastDotPos = position(it);
        ++it;
    }
    if (lastDelimPos == 0)
        return prefix(path, lastDotPos);       
    if (lastDotPos < lastDelimPos)
        return suffix(path, lastDelimPos + 1);
    else
        return infix(path, lastDelimPos + 1, lastDotPos);
}
// ---------------------------------------------------------------------------------------
// Function calculateQPercentile()
// ---------------------------------------------------------------------------------------
// Calculates the p-% percentile of the given values.
unsigned calculatePercentile(String<unsigned> values, const double & p) // No reference on purpose!
{
    SEQAN_ASSERT_LEQ(p, 1.0);
    SEQAN_ASSERT_GT(p, 0.0);
    std::sort(begin(values), end(values));
    unsigned n = length(values);
    unsigned j = floor(n * p);
    bool remainder = j != n * p;
    if(remainder)
        return values[j];
    else
        return values[j-1];
}
// -----------------------------------------------------------------------------
// Function calculatePhredGL()
// -----------------------------------------------------------------------------
// calculate the phred-scaled genotype-likelihood from the given logs of individual genotypel likelihoods.
inline void calculatePhredGL(Call & call, const String<Triple<long double> > & gtLogs)
{
    resize(call.gtLikelihoods, length(gtLogs));
    for (unsigned i = 0; i < length(gtLogs); ++i)
    {
        long double g0 = gtLogs[i].i1;
        long double g1 = gtLogs[i].i2;
        long double g2 = gtLogs[i].i3;
        const long double gTot = log10(exp(g0) + exp(g1) + exp(g2)); // g0-g2 where calculated using ln.
        g0 = -10 * (g0 - gTot);
        g1 = -10 * (g1 - gTot);
        g2 = -10 * (g2 - gTot);
        const long double minPhred = std::min(std::min(g0, g1), g2);
        call.gtLikelihoods[i].i1 = round(g0 - minPhred);
        call.gtLikelihoods[i].i2 = round(g1 - minPhred);
        call.gtLikelihoods[i].i3 = round(g2 - minPhred);
    }
}
// -----------------------------------------------------------------------------
// Function phredsToFrequency()
// -----------------------------------------------------------------------------
// Take the calculatePhredGLs calculated by calculatePhredGL() and estimate the allele frequency of the variant.
inline double phredsToFrequency(const String<Triple<unsigned> > & phredGLs)
{
    unsigned alleles = 0;
    unsigned n = length(phredGLs) * 2;
    for (Iterator<const String<Triple<unsigned> > >::Type it = begin(phredGLs); it != end(phredGLs); ++it)
    {
        if (it->i1 == 0)
        {
            if (it->i2 == 0)        // if the genotype for one sample could not be determined, exclude the sample.
            {
                n -= 2;
            }
        }
        else if (it->i2 == 0)
        {
            alleles += 1;
        }
        else if (it->i3 == 0)
        {
            alleles += 2;
        }
    }
    if (n == 0)
        return 0;
    else
        return static_cast<double>(alleles) / n;
}
#endif /* UTILS_POPDEL_H_ */
