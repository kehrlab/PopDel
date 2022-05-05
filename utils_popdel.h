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
// Class for counting how the tip distances of a read group match the different histograms.
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
    unsigned windowPosition;        // Position of the window.
    unsigned position;              // Assumed start-position (window) of the variant.
    unsigned endPosition;           // Assumed end-position (window) of the variant.
    unsigned significantWindows;   // Number of significant windows that where merged into this variant.
    String<Triple<unsigned> > gtLikelihoods; // PHRED-scaled GT likelihhods of each sample. Order: HomRef, Het, HomDel.
    unsigned char filter;           // 8 Bits indicating the failed filters. 1 at the position indicates failed filter.
    bool hasClipping;               // True if some of the reads leading to the call were soft-clipped.
    bool isOutside;            // True if the window the call is bases on lies outside of the deletion area.
                                    // (right to left) 1.: LR ratio test failed, 2.: High coverage, 3.:Sample Number

    private:
        void init(__uint32 initLen,
                  __uint32 itNum,
                  __uint32 delLen,
                  double llr,
                  double f,
                  unsigned wPos,
                  unsigned pos,
                  unsigned endPos)
        {
            initialLength = initLen;
            iterations = itNum;
            deletionLength = delLen;
            likelihoodRatio = llr;
            frequency = f;
            windowPosition = wPos;
            position = pos;
            endPosition = endPos;
            significantWindows = 0;
            filter = 0;
            hasClipping = false;
            isOutside = false;
        }

    public:
        Call()
        {
            init(0,0,0,0,0,0,0,0);
        }
        Call(__uint32 initLen,
             __uint32 itNum,
             __uint32 delLen,
             double llr,
             double f,
             unsigned wPos,
             unsigned pos,
             unsigned endPos = 0)
        {
            init(initLen, itNum, delLen, llr, f, wPos, pos, endPos);
        }
        Call(__uint32 initLen,
             __uint32 itNum,
             __uint32 delLen,
             double llr,
             double f,
             unsigned wPos,
             Pair<unsigned> startAndEndPos)
        {
            init(initLen, itNum, delLen, llr, f, wPos, startAndEndPos.i1, startAndEndPos.i2);
        }
        void reset(void)
        {
            initialLength = 0;
            iterations = 0;
            deletionLength = 0;
            likelihoodRatio = 0;
            clear(lads);
            clear(dads);
            frequency = 0;
            windowPosition = 0;
            position = 0;
            endPosition = 0;
            significantWindows = 0;
            filter = 0;
            hasClipping = false;
            isOutside = false;
        }
        // =======================================================================================
        // Function print()
        // =======================================================================================
        // Prints some information about the call. Only for debugging and testing purposes,
        void print(void) const
        {
            std::cout << windowPosition << "~" << position << "-" << endPosition << "(" << initialLength << "->" << deletionLength << "%"
            << frequency << "@" << likelihoodRatio << std::endl;
            std::cout << "LADs:";
            for (auto it = begin(lads); it != end(lads); ++it)
                std::cout << *it << ",";
            std::cout << std::endl;
            std::cout << "DADs:";
            for (auto it = begin(dads); it != end(dads); ++it)
                std::cout << "<" << it->ref << "," << it->both << "," << it->between << "," << it->alt
                          << "," << it->right << ">,",
            std::cout << std::endl;
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
// Function setRelWinCovFilter()
// =======================================================================================
// Return true if the realtive window coverage filter has been passed, false otherwise.
inline void setRelWinCovFilter(Call & call)
{
    call.filter |= 16;             // 10000
}
// =======================================================================================
// Function checkRelWinCovPass()
// =======================================================================================
// Return true if the realtive window coverage filter has been passed, false otherwise.
inline bool checkRelWinCovPass(const Call & call)
{
    return !(call.filter & 16);    // 10000
}
// =======================================================================================
// Function checkAllPass()
// =======================================================================================
// Return true if all filters have been passed, false otherwise.
inline bool checkAllPass(const Call & call)
{
    return checkCoveragePass(call) &&
           checkLRPass(call) &&
           checkSamplePass(call) &&
           checkGT0Pass(call) &&
           checkRelWinCovPass(call);
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
inline bool checkDelSizeSimilar(const unsigned & a, const unsigned & b, const double & stddev, const double f = 0.5)
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
    else
    {
        return l >= f * r;
    }
}
inline bool isInDelRange(const Call & a, const Call & b, const double & stddev, const double f = 4)
{
    SEQAN_ASSERT_GEQ(b.position, a.position);
    return (b.position - a.position < (std::min(a.deletionLength, b.deletionLength) + f * stddev));
}
// =======================================================================================
// Function checkAndExtend()
// =======================================================================================
// Return true if the range given by endPos - beginPos of one of the calls needs extension to fit the the estimated del.
// Return false otherwise.
// If true is returned, the endPosition of a is updated to match the endPosition of b.
inline bool checkAndExtend(Call & a, Call & b, const double & stddev)
{
    unsigned aSpan = a.endPosition - a.position;
    unsigned bSpan = b.endPosition - b.position;
    if ((aSpan < a.deletionLength || bSpan < b.deletionLength) && isInDelRange(a, b, stddev))
    {
        a.endPosition = b.endPosition;
        return true;
    }
    else return false;
}
// =======================================================================================
// Function checkEnoughOverlap()
// =======================================================================================
// Return true if the calls a and b overlap for at least for f-percent of the smaller variant's length.
// or if the overlap + 2 stddev is smaller then then smaller variant's length.
inline bool checkEnoughOverlap(const Call & a, const Call & b, const double & stddev, const double f = 0.25)
{
    unsigned aSpan = a.endPosition - a.position;
    unsigned bSpan = b.endPosition - b.position;
    unsigned minLen = std::min(aSpan, bSpan);
    unsigned left = std::max(a.position, b.position);
    unsigned right = std::min(a.position + aSpan, b.position + bSpan);
    int overlap = (right - left);
    if (overlap >= f * minLen)
        return true;
    else if (overlap + 2 * stddev >= minLen)
        return true;
    else
        return false;
}
// =======================================================================================
// Function areSimilar()
// =======================================================================================
// Return true, if the two given SupportStretch are similar enough to be merged.
// If checkAndExtedIs true, the endPosition of a is updated to match the endPosition of b.
inline bool similar(Call & a, Call & b, const double & stddev)
{
    if (checkDelSizeSimilar(a.deletionLength, b.deletionLength, stddev))
    {
        if (checkEnoughOverlap(a, b, stddev))
            return true;
        else
            return checkAndExtend(a, b, stddev);
    }
    else
    {
        return false;
    }
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
// Function setFreqFromGTs()
// =======================================================================================
// Set the allele frequency by counting the number of variant alleles in all samples.
template<typename TCall>
inline void setFreqFromGTs(TCall & call)
{
    unsigned alleleCount = 0;
    for (auto it = begin(call.gtLikelihoods); it != end(call.gtLikelihoods); ++it)
    {
        if (it->i1 == 0)
        {
            continue;
        }
        else if (it->i2 == 0)
        {
            alleleCount += 1;
            continue;
        }
        else
        {
            alleleCount += 2;
        }
    }
    call.frequency = static_cast<double>(alleleCount) / (length(call.gtLikelihoods) * 2);
}
// =======================================================================================
// Function addToLADlist()
// =======================================================================================
// Add the LAD of sample s to the target LAD-list.
inline void addToLADlist(String<Triple<String<unsigned> > > & target,
                         const String<Triple<unsigned> > & source,
                         const unsigned s)
{
    SEQAN_ASSERT_EQ(length(target), length(source));
    SEQAN_ASSERT_LT(s, length(target));
    appendValue(target[s].i1, source[s].i1);
    appendValue(target[s].i2, source[s].i2);
    appendValue(target[s].i3, source[s].i3);
}
// =======================================================================================
// Function addToDADlist()
// =======================================================================================
// Add the LAD of sample s to the target DAD-list.
inline void addToDADlist(String<String<String<unsigned>, Array<5> > > & target,
                         const String<Dad> & source,
                         const unsigned s)
{
    SEQAN_ASSERT_EQ(length(target), length(source));
    SEQAN_ASSERT_LT(s, length(target));
    appendValue(target[s][0], source[s].ref);
    appendValue(target[s][1], source[s].both);
    appendValue(target[s][2], source[s].between);
    appendValue(target[s][3], source[s].alt);
    appendValue(target[s][4], source[s].right);
}
// =======================================================================================
// Function getMedianLAD()
// =======================================================================================
// Return the median LAD and clear the input list of LADs.
inline Triple<unsigned> getMedianLAD(Triple<String<unsigned> > & lads)
{
    unsigned m = length(lads.i1) / 2;
    std::sort(begin(lads.i1), end(lads.i1));
    std::sort(begin(lads.i2), end(lads.i2));
    std::sort(begin(lads.i3), end(lads.i3));
    Triple<unsigned> res = Triple<unsigned>(lads.i1[m], lads.i2[m], lads.i3[m]);
    clear(lads.i1);
    clear(lads.i2);
    clear(lads.i3);
    return res;
}
// =======================================================================================
// Function getMedianDAD()
// =======================================================================================
// Return the median DAD and clear the input list of DADs.
inline Dad getMedianDAD(String<String<unsigned>, Array<5> > & dads)
{
    unsigned m = length(dads[0]) / 2;
    for (unsigned i = 0; i < 5u; ++i)
        std::sort(begin(dads[i]), end(dads[i]));

    Dad res = Dad(dads[0][m], dads[1][m], dads[2][m], dads[3][m], dads[4][m]);
    for (unsigned i = 0; i < 5u; ++i)
        clear(dads[i]);

    return res;
}
// =======================================================================================
// Function skipFailedWindows()
// =======================================================================================
// Advance  currentIt as far as needed to skip over all calls that fail any filters.
// Return false if all calls fail, true otherwise.
inline bool skipFailedCalls(Iterator<String<Call> >::Type & currentIt,
                            const Iterator<const String<Call> >::Type & last)
{
    while (!checkAllPass(*currentIt))        // Skip all windows that did not pass all filters.
    {
        if (currentIt == last)
            return false;
        else
            ++currentIt;
    }
    return currentIt != last;
}
// =======================================================================================
// Function setGenotypes()
// =======================================================================================
// Looks at the windows that are within the deletion length and estimates the genotype for each sample.
// Also sets the LADs and DADs.
// Resets "genotypes"
// Return true on success, false if no genotype could be set (due to lack of windows)
inline bool setGenotypes(const Iterator<String<Call>, Standard >::Type & start,
                         const Iterator<const String<Call>, Standard >::Type end,
                         String<Triple<String<unsigned> > > & lads,
                         String<String<String<unsigned>, Array<5> > > & dads,
                         String<Triple<unsigned> > & genotypes)
{
    unsigned genotypeWinCount = 0;
    for (Iterator<String<Call>, Standard >::Type gtIt = start; gtIt < end; ++gtIt)
    {
        if (gtIt->windowPosition > start->position &&
            gtIt->windowPosition - 30 < start->position + start->deletionLength)
        {   // Only consider genotype likelihoods of windows within the deletion range
            for (unsigned s = 0; s < length(genotypes); ++s)
            {
                genotypes[s].i1 += gtIt->gtLikelihoods[s].i1;
                genotypes[s].i2 += gtIt->gtLikelihoods[s].i2;
                genotypes[s].i3 += gtIt->gtLikelihoods[s].i3;
                addToLADlist(lads, gtIt->lads, s);
                addToDADlist(dads, gtIt->dads, s);
            }
            ++genotypeWinCount;
        }
    }
    if (genotypeWinCount == 0)
        return false;
    // Now get the final genotype likelihoods
    for (unsigned s = 0; s < length(genotypes); ++s)
    {
        double minGt = std::min(std::min(genotypes[s].i1, genotypes[s].i2), genotypes[s].i3);
        double ref = static_cast<double>(genotypes[s].i1 - minGt) / genotypeWinCount;
        double het = static_cast<double>(genotypes[s].i2 - minGt) / genotypeWinCount;
        double hom = static_cast<double>(genotypes[s].i3 - minGt) / genotypeWinCount;
        genotypes[s] = Triple<unsigned> (0, 0, 0);
        start->gtLikelihoods[s].i1 = std::round(ref);
        start->gtLikelihoods[s].i2 = std::round(het);
        start->gtLikelihoods[s].i3 = std::round(hom);
        start->lads[s] = getMedianLAD(lads[s]);
        start->dads[s] = getMedianDAD(dads[s]);
        // Avoid genotypes with nearly equal likelihoods. TODO: A bit hacky. Find better solution.
        if (start->gtLikelihoods[s].i2 == start->gtLikelihoods[s].i3)
        {
            if (het > hom)
                ++(start->gtLikelihoods[s].i2);
            else
                ++(start->gtLikelihoods[s].i3);
        }
        else if (start->gtLikelihoods[s].i1 == start->gtLikelihoods[s].i2)
        {
            if (ref > het)
                ++(start->gtLikelihoods[s].i1);
            else
                ++(start->gtLikelihoods[s].i2);
        }
    }
    return true;
}
// =======================================================================================
// Function getMedianPosition()
// =======================================================================================
// Return the median of the collected position estimates as the final position.
// positions is cleared during this process.
inline unsigned getMedianPosition(String<unsigned> & positions)
{
    SEQAN_ASSERT(!empty(positions));
    std::sort(begin(positions), end(positions));
    unsigned start = positions[length(positions) / 2];
    clear(positions);
    return start;
}
// =======================================================================================
// Function getDeletionLength()
// =======================================================================================
// Return the median of the collected sizes estimates as the final value for the deletion length.
// sizeEstimates is cleared during this process.
inline unsigned getDeletionLength(String<unsigned> & sizeEstimates)
{
    SEQAN_ASSERT(!empty(sizeEstimates));
    std::sort(begin(sizeEstimates), end(sizeEstimates));
    unsigned len = sizeEstimates[length(sizeEstimates) / 2];
    clear(sizeEstimates);
    return len;
}
// =======================================================================================
// Function endPosFallBack()
// =======================================================================================
// Check the clipping of the call and if the endPos matches startPos + deletionLength.
// If there is not sufficient clipping and the positions do not match, fall back to calculating
// EndPos as StarPos + deletionLength
inline void endPosFallBack(unsigned & endPos,
                           const Iterator<String<Call>, Standard >::Type & start,
                           const Iterator<String<Call>, Standard >::Type & last,
                           const unsigned & startPos,
                           const unsigned & deletionLength,
                           const double & stddev)
{
    unsigned lenEnd = startPos + deletionLength;
    if (std::min(lenEnd, endPos) + stddev < std::max(lenEnd, endPos))
    {
        for (Iterator<String<Call>, Standard >::Type it = start; it  < last; ++it)
        {
            if (it->hasClipping)
                return;
        }
        //std::cout << "Falling Back: " << startPos << "-" << endPos << " (" << deletionLength << ") " << "\t->\t" << lenEnd << " (DIFF: " << static_cast<int>(endPos) - static_cast<int>(lenEnd) << ")" <<std::endl;
        endPos = lenEnd;
    }
}
inline void mergeWindowRange(const Iterator<String<Call>, Standard >::Type & start,
                             const Iterator<String<Call>, Standard >::Type & last,
                             String<Triple<String<unsigned> > > & lads,
                             String<String<String<unsigned>, Array<5> > > & dads,
                             String<Triple<unsigned> > & genotypes,
                             String<unsigned> & startPositions,
                             String<unsigned> & endPositions,
                             String<unsigned> & sizeEstimates,
                             long double & lr,
                             unsigned & callCount,
                             unsigned & winCount,
                             unsigned & significantWindows,
                             const double & r,
                             const double & stddev)
{
    start->position = getMedianPosition(startPositions);
    start->endPosition = getMedianPosition(endPositions);
    start->deletionLength = getDeletionLength(sizeEstimates);
    endPosFallBack(start->endPosition, start, last, start->position, start->deletionLength, stddev);
    start->likelihoodRatio = lr / winCount;
    if (setGenotypes(start, last, lads, dads, genotypes))
    {
        setFreqFromGTs(*start);

        start->significantWindows = significantWindows;
        if (30.0 * significantWindows / start->deletionLength < r)
            setRelWinCovFilter(*start);

        winCount = 1;
        significantWindows = 1;
        lr = 0.0;
        ++callCount;
    }
    else
    {
        markInvalidCall(*start);
    }
}
// =======================================================================================
// Function mergeDeletionWindows()
// =======================================================================================
// Merge window deletion calls.
// Return false if the string of calls is empty, true otherwise.
inline bool mergeDeletionWindows(String<Call> & calls,
                                 const double & stddev,
                                 const double & r,
                                 const bool outputFailed = false)
{
    if (length(calls) <= 1u)
        return false;

    std::sort(begin(calls), end(calls), lowerCall);
    Iterator<String<Call> >::Type currentIt = begin(calls, Standard());
    const Iterator<String<Call>, Standard >::Type last = end(calls) - 1;
    if (!outputFailed)
        if (!skipFailedCalls(currentIt, last))
            return false;

    const Iterator<String<Call> >::Type firstGoodWin = currentIt;
    String<Triple<unsigned> > genotypes;
    resize(genotypes, length(currentIt->gtLikelihoods), Triple<unsigned>(0, 0, 0));
    String<unsigned> startPositions;
    String<unsigned> endPositions;
    String<unsigned> sizeEstimates;
    String<Triple<String<unsigned> > > lads;
    resize(lads, length(currentIt->lads));
    String<String<String<unsigned>, Array<5> > > dads;
    resize(dads, length(currentIt->dads));
    for (Iterator<String<String<String<unsigned>, Array<5> > > >::Type it = begin(dads); it != end(dads); ++it)
        resize(*it, 5, Exact());

    append(startPositions, currentIt->position);
    append(endPositions, currentIt->endPosition);
    append(sizeEstimates, currentIt->deletionLength);
    unsigned callCount = 1;
    unsigned winCount = 1;
    unsigned significantWindows = 1;
    long double lr = currentIt->likelihoodRatio;

    Iterator<String<Call>, Standard >::Type it = firstGoodWin + 1;
    while (true)
    {
        if (similar(*currentIt, *it, stddev))
        {
            if (checkAllPass(*it))
            {
                appendValue(startPositions, it->position);
                append(endPositions, it->endPosition);
                appendValue(sizeEstimates, it->deletionLength);
                ++significantWindows;
            }
            ++winCount;
            lr += it->likelihoodRatio;
            markInvalidCall(*it);

            if (it == last)
            {
                if (!empty(startPositions))
                    mergeWindowRange(currentIt,
                                     it,
                                     lads,
                                     dads,
                                     genotypes,
                                     startPositions,
                                     endPositions,
                                     sizeEstimates,
                                     lr,
                                     callCount,
                                     winCount,
                                     significantWindows,
                                     r,
                                     stddev);

                break;
            }
        }
        else
        {
            if (winCount != 1 && !empty(startPositions))
                mergeWindowRange(currentIt,
                                 it,
                                 lads,
                                 dads,
                                 genotypes,
                                 startPositions,
                                 endPositions,
                                 sizeEstimates,
                                 lr,
                                 callCount,
                                 winCount,
                                 significantWindows,
                                 r,
                                 stddev);
            else
                markInvalidCall(*currentIt);

            currentIt = it;
        }
        if (it != last)
        {
            ++it;
        }
        else
        {
            if (winCount == 1)
            {
                --callCount;
                markInvalidCall(*currentIt);
            }
            break;
        }
    }
    String<Call> tmp;
    reserve(tmp, callCount, Exact());
    it = firstGoodWin;
    for (; it <= last; ++it)
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
inline double sum(const Triple<double> & t)
{
    return (t.i1 + t.i2 + t.i3);
}
inline long double sum(const Triple<long double> & t)
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
// Function setNoDataSamples()
// =======================================================================================
// Sets all genotypes for samples marked in lowCoverageSamples to 0.
inline void setNoDataSamples(Call & currentCall, const String<bool> & lowCoverageSamples)
{
    for (unsigned i = 0; i < length(lowCoverageSamples); ++i)
    {
        if (lowCoverageSamples[i])
        {
            currentCall.gtLikelihoods[i].i1 = 0;
            currentCall.gtLikelihoods[i].i2 = 0;
            currentCall.gtLikelihoods[i].i3 = 0;
        }
    }
}
// =======================================================================================
// Function checkSampleNumber()
// =======================================================================================
// Checks if the number of samples with data at the current position is high enough.
// Return true if so, false otherwise.
inline bool checkSampleNumber(const String<bool> & lowCoverageSamples, const double & threshold)
{
    SEQAN_ASSERT(!empty(lowCoverageSamples));
    unsigned total = length(lowCoverageSamples);
    unsigned dataSample = total;
    for (unsigned  i = 0; i < total; ++i)
        if (lowCoverageSamples[i])
            --dataSample;

    return static_cast<double>(dataSample) / total >= threshold;
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

template<typename TNum>
inline TNum max(const Triple<TNum> & t)
{
    TNum currentMax = t.i1 >= t.i2 ? t.i1:t.i2;
    if (t.i3 > currentMax)
        currentMax = t.i3;

    return currentMax;
}
template<typename TNum>
inline unsigned whichMax(const Triple<TNum> & t)
{
    if (t.i1 >= t.i2)
    {
        if (t.i1 >= t.i3)
            return 0;
        else
            return 2;
    }
    else if (t.i2 >= t.i3)
    {
        return 1;
    }
    else
        return 2;
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
// Function getIndexFileName()
// =======================================================================================
// Check if the alignment file is a BAM or CRAM file and append .bai or .crai accoordingly.
inline String<char> getIndexFileName(const HtsFile & file)
{
    const htsFormat & format = file.fp->format;
    String<char> indexFileName = file.filename;
    if (format.format == htsExactFormat::bam)
    {
        append(indexFileName, ".bai");
    }
    else if (format.format == htsExactFormat::cram)
    {
        append(indexFileName, ".crai");
    }
    else
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not determine file type from file name for \'" << file.filename << "\'."
            <<  " Please make sure that the file name ends with \'.bam\' or \'.cram\'. Terminating.";
        SEQAN_THROW(IOError(toCString(msg.str())));
    }
    return indexFileName;
}

// =======================================================================================
// Function loadBai()
// =======================================================================================
// Load the index belonging to the given BAM/CRAM-file
// Return false on errors, true otherwise.
inline bool loadBaiCrai(HtsFileIn & infile)
{
    String<char> indexPathIn = getIndexFileName(infile);
    if (!loadIndex(infile, toCString(indexPathIn)))
    {
        CharString message = "[PopDel] ERROR: Could not load the index file \'";
        message += indexPathIn;
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
    reg = GenomicRegion();   // Necessary to reset reg, because parse() does not guarantee to overwrite all old entries.
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
// Function getReadGroup()
// =======================================================================================
// Extract the read-group encoded in tags.
// Return the read-group as a CharString.
inline CharString getReadGroup(const CharString & tags)
{
    BamTagsDict dict(tags);                                        // TODO: Catch if RG tag doen't exist. Error message.
    unsigned key = 0;
    findTagKey(key, dict, "RG");                                   //TODO faster access if position in dict is known
    CharString rg = "";
    extractTagValue(rg, dict, key);
    return rg;
}
// =======================================================================================
// Function getReadGroup()
// =======================================================================================
// Extract the read-group encoded in tags
// Return its rank in the header by extracting it from in the map of read groups.
inline unsigned getReadGroup(const CharString & tags,
                             const std::map<CharString, unsigned> & readGroups,
                             bool merge = false)
{
    if (merge)
        return readGroups.begin()->second;

    CharString rg = getReadGroup(tags);
    SEQAN_ASSERT_NEQ(readGroups.count(rg), 0u);
    return readGroups.at(rg);
}
// =======================================================================================
// Function getRgIdFromKstring()
// =======================================================================================
// Extact the ID from the kstring of the @RG line in the header.
// Return true on success, false otherwise.
inline bool getRgIdFromKstring(CharString & id, const kstring_t & k)
{
    SEQAN_ASSERT(empty(id));
    SEQAN_ASSERT_GT(k.l, 7u);
    SEQAN_ASSERT_EQ(k.s[0], '@');
    SEQAN_ASSERT_EQ(k.s[1], 'R');
    SEQAN_ASSERT_EQ(k.s[2], 'G');
    SEQAN_ASSERT_EQ(k.s[3], '\t');
    SEQAN_ASSERT_EQ(k.s[4], 'I');
    SEQAN_ASSERT_EQ(k.s[5], 'D');
    SEQAN_ASSERT_EQ(k.s[6], ':');
    unsigned i = 7;
    while (i < k.l && k.s[i] != '\t')
    {
        appendValue(id, k.s[i]);
        ++i;
    }
    return !empty(id);
}
// =======================================================================================
// Function getRgIdFromKstring()
// =======================================================================================
// Extact the SM from the kstring of the @RG line in the header.
// Return true on success, false otherwise.
inline bool getSMFromKstring(CharString & sm, const kstring_t & k)
{
    SEQAN_ASSERT(empty(sm));
    SEQAN_ASSERT_GT(k.l, 7u);
    SEQAN_ASSERT_EQ(k.s[0], '@');
    SEQAN_ASSERT_EQ(k.s[1], 'R');
    SEQAN_ASSERT_EQ(k.s[2], 'G');
    SEQAN_ASSERT_EQ(k.s[3], '\t');
    SEQAN_ASSERT_EQ(k.s[4], 'I');
    SEQAN_ASSERT_EQ(k.s[5], 'D');
    SEQAN_ASSERT_EQ(k.s[6], ':');
    unsigned i = 9; // 7 Is the first letter of the ID, 8 or above is the tab after the ID
    clear(sm);
    while (i < k.l && empty(sm))
    {
        if (k.s[i] == 'S' && i < k.l - 1)
        {
            ++i;
            if (k.s[i] == 'M' && i < k.l - 1)
            {
                ++i;
                if (k.s[i] == ':' && i < k.l - 1)
                {
                    ++i;
                    while (i < k.l && k.s[i] != '\t')
                    {
                        appendValue(sm, k.s[i]);
                        ++i;
                    }
                }
            }
        }
        ++i;
    }
    return !empty(sm);
}
// =======================================================================================
// Function getReadGroups()
// =======================================================================================
// Extract all read group IDs and their rank in the header and write them to map.
// Return the number of read groups in the header.
inline unsigned getReadGroups(std::map<CharString, unsigned> & readGroups, HtsFile & file, bool merge = false)
{
    readGroups.clear();
    int good = 0;
    int pos = 0;
    while (true)
    {
        kstring_t kstring = KS_INITIALIZE;
        good = sam_hdr_find_line_pos(file.hdr, "RG", pos, &kstring);
        if (good == 0)
        {
            CharString id = "";
            if (getRgIdFromKstring(id, kstring))
            {
                if (merge)
                    readGroups[id] = 0u;
                else
                    readGroups[id] = static_cast<unsigned>(pos);

                ++pos;
            }
            else
            {
                std::ostringstream msg;
                msg << "[PopDel] Error while trying to parse ID from @RG tag. Terminating";
                SEQAN_THROW(ParseError(toCString(msg.str())));
            }
        }
        else
        {
            ks_free(&kstring);
            break;
        }
        ks_free(&kstring);
    }
    return readGroups.size();
}
// =======================================================================================
// Function getSampleName()
// =======================================================================================
// Get the sample name encoded in the @RG SM-tag.
// Return true on success, false otherwise.
// Throw an error, if multiple conflicting sample names are found.
inline bool getSampleName(CharString & sampleName, HtsFile & file)
{
    int good = 0;
    int pos = 0;
    while (true)
    {
        kstring_t kstring = KS_INITIALIZE;
        good = sam_hdr_find_line_pos(file.hdr, "RG", pos, &kstring);
        if (good == 0)
        {
            CharString sm = "";
            if (getSMFromKstring(sm, kstring))
            {
                if (sampleName == "")
                {
                    sampleName = sm;
                }
                else if (sampleName != sm)
                {
                    std::ostringstream msg;
                    msg << "[PopDel] Found conflicting SM-tags in the alignment file! Previsously found 'SM:"
                        << sampleName
                        << "'. But now found 'SM:"
                        << sm
                        << "'. PopDel profile must only be used on single sample alignemnt files. Terminating.";
                    SEQAN_THROW(ParseError(toCString(msg.str())));
                }
            }
            ++pos;
        }
        else
        {
            ks_free(&kstring);
            break;
        }
        ks_free(&kstring);
    }
    return !empty(sampleName);
}
// =======================================================================================
// Functions for genomic interval loading and processing
// =======================================================================================
// =======================================================================================
// Function _getWholeGenomeIntervals()
// =======================================================================================
// Append start and end positions + sequence names of all contigs/chromosomes and store them in String of GenomicRegion.
// Return total length of all intervals
inline unsigned getWholeGenomeIntervals(String<GenomicRegion> & intervals, const HtsFile & infile)
{
    String<String<char> > contigNames;
    String<int> contigLengths;
    getContigNames(contigNames, infile);
    getContigLengths(contigLengths, infile);
    SEQAN_ASSERT_EQ(length(contigNames), length(contigLengths));
    unsigned totalLength = 0;
    auto lenIt = begin(contigLengths);
    for (auto it = begin(contigNames); it != end(contigNames); ++it, ++lenIt)
    {
        GenomicRegion itv;
        itv.seqName = *it;                                          //Store chromosome/contig name of region
        itv.beginPos = 0;
        itv.endPos = *lenIt;                                        //Get the end position of the interval
        appendValue(intervals, itv);                                //Append new interval to list of intervals
        totalLength += itv.endPos;                                  //Increase total length of sum of intervals
    }
    return totalLength;
}
// =======================================================================================
// Struct lowerGenomicRegion()
// =======================================================================================
// Return true if the first GenomicRegion starts before the second one starts.
// Note that this function is based on the lexicographical order of the seq names and therefore must NOT be used
// when the comparison is anyhow related to the order of the contigs in the bam file or profile!.
// Meant for stable_sort.
inline bool lowerGenomicRegion(const GenomicRegion & r1, const GenomicRegion & r2)
{
    Lexical<> cmp(r1.seqName, r2.seqName);                      //Lexicographically compare the sequence names.
    if (isLess(cmp))                                            //r1.seqName is lexicographically smaller
        return true;
    else if (isGreater(cmp))                                    //r1.seqName is lexicographically bigger
        return false;
    else
        return r1.beginPos < r2.beginPos;
}
// =======================================================================================
// Struct lowerRIDGenomicRegion()
// =======================================================================================
// Return true if the first GenomicRegion starts before the second one starts.
// This Function uses the rIDs for defining the ordering of contigs.
inline bool lowerRIDGenomicRegion(const GenomicRegion & r1, const GenomicRegion & r2)
{
    if (r1.rID < r2.rID)                                            //r1.seqName is lexicographically smaller
        return true;
    else if (r1.rID > r2.rID)                                    //r1.seqName is lexicographically bigger
        return false;
    else
        return r1.beginPos < r2.beginPos;
}
// =======================================================================================
// Function _fillInvalidPositions()
// =======================================================================================
// Replace invalid values of GenomicRegions object with processable values.
// Use 0 for start position. Use end position of sequence for end position of interval.
inline void fillInvalidPositions(GenomicRegion & itv,
                                 const String<String<char> > & contigNames,
                                 const String<int32_t> & contigLengths)
{
    if (itv.beginPos == GenomicRegion::INVALID_POS)
        itv.beginPos = 0;
    if (itv.endPos == GenomicRegion::INVALID_POS)
    {
        for (unsigned i = 0; i < length(contigNames); ++i)
        {
            if (contigNames[i] == itv.seqName)
                itv.endPos = contigLengths[i];
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
void createRegularIntervals(String<GenomicRegion> & intervals,
                            const String<GenomicRegion> & rois,
                            const CharString & referenceVersion)
{
    String<GenomicRegion> tmpIntervals;
    resize(tmpIntervals, 22, Exact());
    std::ostringstream msg;
    CharString r = referenceVersion;
    toLower(r);
    std::string prefix = "";
    bool oldBuild = true;

    if (r == "grch38" || r == "hg38")
    {
        prefix = "chr";
        oldBuild = false;
    }
    else if ( r == "hg19")
    {
        prefix = "chr";
    }
    else if  (r == "grch37")
    {
        prefix = "";
    }

    if (!oldBuild)
    {
        msg << "Using default sampling intervals for GRCh38/hg38.";
        parse(tmpIntervals[0], prefix + "1:35000000-36000000");
        parse(tmpIntervals[1], prefix + "2:174000000-175000000");
        parse(tmpIntervals[2], prefix + "3:36500000-37500000");
        parse(tmpIntervals[3], prefix + "4:88000000-89000000");
        parse(tmpIntervals[4], prefix + "5:38000000-39000000");
        parse(tmpIntervals[5], prefix + "6:38000000-39000000");
        parse(tmpIntervals[6], prefix + "7:38000000-39000000");
        parse(tmpIntervals[7], prefix + "8:19000000-20000000");
        parse(tmpIntervals[8], prefix + "9:19000000-20000000");
        parse(tmpIntervals[9], prefix + "10:19000000-20000000");
        parse(tmpIntervals[10], prefix + "11:19000000-20000000");
        parse(tmpIntervals[11], prefix + "12:19000000-20000000");
        parse(tmpIntervals[12], prefix + "13:25000000-26000000");
        parse(tmpIntervals[13], prefix + "14:25000000-26000000");
        parse(tmpIntervals[14], prefix + "15:25000000-26000000");
        parse(tmpIntervals[15], prefix + "16:25000000-26000000");
        parse(tmpIntervals[16], prefix + "17:31000000-32000000");
        parse(tmpIntervals[17], prefix + "18:31000000-32000000");
        parse(tmpIntervals[18], prefix + "19:31000000-32000000");
        parse(tmpIntervals[19], prefix + "20:33000000-34000000");
        parse(tmpIntervals[20], prefix + "21:21000000-22000000");
        parse(tmpIntervals[21], prefix + "22:25000000-26000000");
    }
    else
    {
        msg << "Using default sampling intervals for GRCh37/hg19.";
        parse(tmpIntervals[0], prefix + "1:35000000-36000000");
        parse(tmpIntervals[1], prefix + "2:174000000-175000000");
        parse(tmpIntervals[2], prefix + "3:36500000-37500000");
        parse(tmpIntervals[3], prefix + "4:88000000-89000000");
        parse(tmpIntervals[4], prefix + "5:38000000-39000000");
        parse(tmpIntervals[5], prefix + "6:38000000-39000000");
        parse(tmpIntervals[6], prefix + "7:37000000-38000000");
        parse(tmpIntervals[7], prefix + "8:19000000-20000000");
        parse(tmpIntervals[8], prefix + "9:19000000-20000000");
        parse(tmpIntervals[9], prefix + "10:19000000-20000000");
        parse(tmpIntervals[10], prefix + "11:19000000-20000000");
        parse(tmpIntervals[11], prefix + "12:19000000-20000000");
        parse(tmpIntervals[12], prefix + "13:30000000-31000000");
        parse(tmpIntervals[13], prefix + "14:26000000-27000000");
        parse(tmpIntervals[14], prefix + "15:26000000-27000000");
        parse(tmpIntervals[15], prefix + "16:25000000-26000000");
        parse(tmpIntervals[16], prefix + "17:31000000-32000000");
        parse(tmpIntervals[17], prefix + "18:31000000-32000000");
        parse(tmpIntervals[18], prefix + "19:31000000-32000000");
        parse(tmpIntervals[19], prefix + "20:33000000-34000000");
        parse(tmpIntervals[20], prefix + "21:21000000-22000000");
        parse(tmpIntervals[21], prefix + "22:34000000-35000000");
    }
    printStatus(msg);
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
               " the parameter estimation. If you are using hg19/GRCh37 please use '--reference', followed by your reference build."
               " Please note, that your contig names have to match the patterns definde by used by the builds."
               " Otherwise, please check the contig names of the ROI's and/or use user-defined"
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
                          const HtsFile & infile,
                          const String<GenomicRegion> & rois,
                          const CharString & referenceVersion)
{
    if (filename != "")
    {
        size_t initialLength = length(intervals);
        std::ifstream intervalFile(toCString(filename));
        if (!intervalFile.is_open())
        {
            std::ostringstream msg;
            msg << "[PopDel] Could not open interval file \'" << filename << "\' for reading.";
            SEQAN_THROW(IOError(toCString(msg.str())));
        }
        String<String<char> > contigNames;
        String<int> contigLengths;
        getContigNames(contigNames, infile);
        getContigLengths(contigLengths, infile);
        std::string word;
        GenomicRegion itv;
        while (intervalFile >> word)
        {
            parseGenomicRegion(itv, word);                                                           //rIDs will not be set
            fillInvalidPositions(itv, contigNames, contigLengths);
            appendValue(intervals, itv);
        }
        std::ostringstream msg;
        msg << "Finished reading " << (length(intervals) - initialLength) << " intervals from file \'" << filename << "\'.";
        printStatus(msg);
    }
    else
    {
        createRegularIntervals(intervals, rois, referenceVersion);
    }
}
// =======================================================================================
// Function expandIntervals()
// =======================================================================================
// Expands the intervals by the deletionSize.
inline void expandIntervals(String<GenomicRegion> & intervals, const int & maxDeletionLength)
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
    HtsFileIn infile;
    open(infile, toCString(bamfile));
    String<CharString> contigNames;
    String<int> contigLengths;
    getContigNames(contigNames, infile);
    getContigLengths(contigLengths, infile);
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
        fillInvalidPositions(itv, contigNames, contigLengths);
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
    reserve(mergedIntervals, length(intervals), Exact());
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
inline void setRIDs(String<GenomicRegion> & intervals, const HtsFile & infile)
{                                                   //TODO: Check if regions without valid ID is realy skipped lateron
    typedef Iterator<String<GenomicRegion> >::Type TItvIter;
    std::map<CharString, int32_t> contigNameMap;
    getContigNameToIDMap(contigNameMap, infile);
    SEQAN_ASSERT_GT(contigNameMap.size(), 0u);
    TItvIter itv = begin(intervals);
    TItvIter itvEnd = end(intervals);
    while (itv != itvEnd)
    {
        auto it = contigNameMap.find(toCString(itv->seqName));
        if (it != contigNameMap.end())
        {
            itv->rID = it->second;
        }
        else
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
// Function rIDAlreadyProcessed()
// =======================================================================================
// Return true if the rID has been marked as completely processed in finishedRIDs.
// Return false otherwise.
inline bool rIDAlreadyProcessed(const String<bool> & finishedRIDs, const int32_t rID)
{
    SEQAN_ASSERT_LT(rID, static_cast<int32_t>(length(finishedRIDs)));
    return finishedRIDs[rID];
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
// Function getSampleNameFromPath()
// =======================================================================================
// Extract the last part of the path to get the file/sample name without file format ending.
inline Infix<const String<char> >::Type getSampleNameFromPath(const String<char> & path)
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
struct Sorted{};
// -----------------------------------------------------------------------------
// Function getPercentile()
// -----------------------------------------------------------------------------
// Return the value in s at percentile p. s will be sorted in the process, if the flag Sorted() ist not useed
// Uses the nearest rank method: No more then p*100 percent of the values in s are < returned value.
// Note: Don't try to use this on an (potentially) empty string!
inline unsigned getPercentile(const String<unsigned> & s, double p, Sorted sorted)
{
    (void) sorted;
    SEQAN_ASSERT_GT(length(s), 0u);
    int lPre = std::ceil(static_cast<double>(length(s)) * p);
    unsigned l = lPre == 0 ? 0u : lPre - 1;
    return s[l];
}
inline unsigned getPercentile(String<unsigned> & s, double p = 0.8)
{
    std::sort(begin(s), end(s));
    return getPercentile(s, p, Sorted());
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
inline bool addToSampleNames(String<CharString> & sampleNames, const CharString & sampleName)
{
    CharString s = sampleName;
    unsigned c = 1;
    for (Iterator<const String<CharString> >::Type it = begin(sampleNames); it != end(sampleNames); ++it)
    {
        if (*it == s)
        {
            ++c;
            s = sampleName;
            appendValue(s, '#');
            append(s, std::to_string(c));
        }
    }
    if (c > 1)
    {
        std::ostringstream msg;
        msg << "WARNING: Duplicate sample name '" << sampleName << "'. Renaming it to " << s << ".";
        printStatus(msg);
    }
    appendValue(sampleNames, s);
    return c == 1;
}
// Return a Pair consisting of the left and right border that determines if a read supports a deletion or not.
inline Pair<int> getDelSupportBorders(const unsigned & deletion_length,
                                      const int & lowerQuantileDist,
                                      const int & upperQuantileDist)
{
    return Pair<int>(deletion_length - lowerQuantileDist, deletion_length + upperQuantileDist);
}
// Return true if the given insert size deviatino supports the given deletion, false otherwise.
inline bool supportsDel(const int & deviation, const Pair<int> & borders)
{
    return deviation >= borders.i1 && deviation <= borders.i2;
}
// -----------------------------------------------------------------------------
// Function add()
// -----------------------------------------------------------------------------
// Add the values in r to those in l.
template <typename TNumL, typename TNumR>
inline void add(String<TNumL> & l, const String<TNumR> & r)
{
    SEQAN_ASSERT_EQ(length(l), length(r));
    if (empty(r))
        return;

    for (unsigned i = 0; i < length(l); ++i)
        l[i] += r[i];
}
#endif /* UTILS_POPDEL_H_ */
