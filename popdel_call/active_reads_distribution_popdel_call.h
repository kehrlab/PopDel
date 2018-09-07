#ifndef DISTRIBUTION_MODEL_POPDEL_H_
#define DISTRIBUTION_MODEL_POPDEL_H_

#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "profile_structure_popdel_call.h"
#include "parameter_parsing_popdel_call.h"

using namespace seqan;

// ======================================================================================
// Class SimulationParameters
// ======================================================================================
//Struct which holds all parameters necessary for the simulation of active read pairs around deletions.
struct SimulationParameters
{
    String<double> & cumulativeDensity;
    unsigned windowSize;
    unsigned deletionSize;
    unsigned readSize;
    double avgNewActiveReadsPerWin;
    unsigned delStartWin;
    double alleleFrequency;
    unsigned offset;
    unsigned flankingWindows;

    SimulationParameters(String<double>& ref,
                         unsigned w,
                         unsigned d, 
                         unsigned r,
                         double a,
                         double h,
                         unsigned o,
                         unsigned f):
    cumulativeDensity(ref),
    windowSize(w),
    deletionSize(d),
    readSize(r),
    avgNewActiveReadsPerWin(a),
    delStartWin(f + 1),
    alleleFrequency(h),
    offset(o),
    flankingWindows(f){}
};
// ======================================================================================
// Tags
// ======================================================================================
// Tag for indicating that the expected reads should be added to the given String.
struct Add
{};
struct FirstHalf
{};
// ======================================================================================
// Class SupportCandidate
// ======================================================================================
struct SupportCandidate
{
    __uint32 start;
    __uint32 end;                   // One window behind the last window.
    __uint32 delSize;
    double   frequency;
    double   bestLikelihoodRatio;   // best likelihood ratio of a single window in the candidate
    bool     fresh;                 // true, if start may be set anew.

    SupportCandidate():
    start(maxValue<__uint32>()),
    end(maxValue<__uint32>()),
    delSize(0),
    frequency(0.0),
    bestLikelihoodRatio(0),
    fresh(true){};
};
// ======================================================================================
// Class SupportStretch
// ======================================================================================
struct SupportStretch
{
    __uint32 start;
    __uint32 numWindows;          // Number of windows the stretch extendes (inlcuding the first one)
    __uint32 delSize;
    __uint32 bestWindow;          // Most likely starting window of the deletion.
    double   frequency;
    double   bestLikelihoodRatio;

    SupportStretch():
    start(maxValue<__uint32>()),
    numWindows(maxValue<__uint32>()),
    delSize(0),
    bestWindow(maxValue<__uint32>()),
    frequency(0),
    bestLikelihoodRatio(0)
    {};
    SupportStretch(unsigned s, unsigned e, unsigned d = 0, unsigned b = maxValue<__uint32>(), double f = 0):
    start(s), numWindows(e), delSize(d), bestWindow(b), frequency(f)
    {};
    SupportStretch(const SupportCandidate & c, const unsigned & windowSize):
    start(c.start),
    numWindows((c.end - c.start) / windowSize),   //c.end is one window behind the last window with signal/valid length.
    delSize(c.delSize),
    bestWindow(maxValue<__uint32>()),
    frequency(c.frequency),
    bestLikelihoodRatio(c.bestLikelihoodRatio)
    {};
};
// ---------------------------------------------------------------------------------------
// Function checkCombinedConsistency()
// ---------------------------------------------------------------------------------------
// Check if the consistency (i.e. the ratio of numWindows and deletioSize) of the combined variants is higher than
// the mean consistency of the individual variants.
inline bool checkCombinedConsistency(const SupportStretch & a, const SupportStretch & b, unsigned windowSize)
{
    double aDelWin = static_cast<double>(a.delSize) / windowSize;
    double bDelWin = static_cast<double>(b.delSize) / windowSize;
    double aConsistency = std::min(aDelWin, static_cast<double>(a.numWindows)) /
                          std::max(aDelWin, static_cast<double>(a.numWindows));
    double bConsistency = std::min(bDelWin, static_cast<double>(b.numWindows)) /
                          std::max(bDelWin, static_cast<double>(b.numWindows));
    unsigned leftStart= std::min(a.bestWindow, b.bestWindow);
    unsigned rightEnd= std::max(a.bestWindow + a.numWindows * windowSize, b.bestWindow + b.numWindows * windowSize);
    double combinedWinSize = rightEnd - leftStart;
    double meanDelSize = static_cast<double>(a.delSize + b.delSize) / 2; // TODO: Maybe weight by LR
    double combinedConsistency = std::min(combinedWinSize , meanDelSize) / std::max(combinedWinSize, meanDelSize);
    return (combinedConsistency > (aConsistency + bConsistency) / 2 + 0.0000000001);// + e-10 to counter precision loss.
}
// ---------------------------------------------------------------------------------------
// Function checkEnoughOverlap()
// ---------------------------------------------------------------------------------------
// Return true if the variants a and b overlap for at least for f-percent of the smaller variant's number of windows.
inline bool checkEnoughOverlap(const SupportStretch & a, const SupportStretch & b, const unsigned & windowSize, const double f = 0.5)
{
    unsigned minLen = std::min(a.numWindows, b.numWindows) * windowSize;
    unsigned left = std::max(a.bestWindow, b.bestWindow);
    unsigned right = std::min(a.bestWindow + a.numWindows * windowSize, b.bestWindow + b.numWindows * windowSize);
    int overlap = (right - left);
    return overlap >= f * minLen;
}
// ---------------------------------------------------------------------------------------
// Function areSimilar()
// ---------------------------------------------------------------------------------------
// Return true, if the two given SupportStretch are similar enough to be considered equal.
inline bool areSimilar(const SupportStretch & a, const SupportStretch & b, const unsigned & windowSize, const double & stddev)
{
    if (checkCombinedConsistency(a, b, windowSize))
        return true;
    else
        return checkDelSizeSimilar(a.delSize, b.delSize, stddev) && checkEnoughOverlap(a, b, windowSize);
}
// ---------------------------------------------------------------------------------------
// Function mergeStretches()
// ---------------------------------------------------------------------------------------
// Merge 2 supportStretches into the first one.
inline void mergeStretches(SupportStretch & a, SupportStretch & b, const unsigned & windowSize)
{
    a.start = std::min(a.start, b.start);
    if (a.bestLikelihoodRatio < b.bestLikelihoodRatio)
    {
        a.delSize = b.delSize;
        a.bestLikelihoodRatio = b.bestLikelihoodRatio;
    }
    a.numWindows = std::ceil(static_cast<double>(a.delSize) / windowSize);  // TODO: Examine how this behaves.
}
// ======================================================================================
// Class SupportChangeMap
// ======================================================================================
struct SupportChangeMap
{
    String<SupportStretch> supportStretches;                   // Has to be ordered!

    SupportChangeMap(): supportStretches()
    {
        reserve(supportStretches, 1000);
    };
    // ======================================================================================
    // Function add()
    // ======================================================================================
    inline void add(const unsigned & s, const unsigned & e)
    {
        appendValue(supportStretches, SupportStretch(s, e));
    }
    inline void add(const SupportStretch & s)
    {
        appendValue(supportStretches, s);
    }
    inline void add(SupportCandidate & c, const unsigned & windowSize)
    {
        if (c.end != maxValue<__uint32>() && (c.end - c.start) / windowSize >= 2u)
        {
            appendValue(supportStretches, SupportStretch(c, windowSize));
            --(back(supportStretches).numWindows);                          //Close the window one position earlier.
            c.start = maxValue<__uint32>();
            c.end = maxValue<__uint32>();
            c.delSize = 0;
            c.frequency = 0;
        }
    }
    // ======================================================================================
    // Function merge()
    // ======================================================================================
    // TODO: Write Test. // TODO: Implement mor efficient alternative.
    inline void merge(SupportChangeMap & m, const unsigned & windowSize, const double & stddev)
    {
        Iterator<String<SupportStretch>, Rooted>::Type itA = begin(supportStretches, Rooted());
        Iterator<String<SupportStretch>, Rooted>::Type itB = begin(m.supportStretches, Rooted());
        String<SupportStretch> tmp;
        reserve(tmp, length(supportStretches) + length(m.supportStretches));
        while(!atEnd(itA) && !atEnd(itB))
        {
            bool merge = areSimilar(*itA, *itB, windowSize, stddev);
            if (merge)
            {
                mergeStretches(*itA, *itB, windowSize);
                appendValue(tmp, *itA);
                ++itA;
                ++itB;
            }
            else
            {
                if (itA->start < itB->start)                              // A starts before B
                {
                    appendValue(tmp, *itA);
                    ++itA;
                }
                else if (itA->start == itB->start)                        // They have the same start
                {
                    if (itA->numWindows < itB->numWindows)                // but A ends before B
                    {
                        appendValue(tmp, *itA);
                        ++itA;
                    }
                    else if (itA->numWindows == itB->numWindows) // They are the same. This should not happen anymore.
                    {
                        mergeStretches(*itA, *itB, windowSize);
                        ++itA;
                        ++itB;
                    }
                    else                                                   // Same start, but B ends before A
                    {
                        appendValue(tmp, *itB);
                        ++itB;
                    }
                }
                else                                                       // B starts before A
                {
                    appendValue(tmp, *itB);
                    ++itB;
                }
            }
        }
        //Append the rest.
        while(!atEnd(itA))
        {
            appendValue(tmp, *itA);
            ++itA;
        }
        while(!atEnd(itB))
        {
            appendValue(tmp, *itB);
            ++itB;
        }
        move(supportStretches, tmp);
        resize(m.supportStretches, 0);
    }
};
// ======================================================================================
// Function stretchTooLong()
// ======================================================================================
// Return true if the distance between pos and start is too big.
inline bool stretchTooLong(const unsigned & pos, const __uint32 & start, const unsigned & delSize)
{
    return pos - start > 2 * delSize;
}
// ======================================================================================
// Function calculateLimits()
// ======================================================================================
// Return minDev and maxDev (as deviations from mean).
// Remember: hist.lowerQuantile and hist.upperQuantile are the distances from the mean to the qantile borders.
inline Pair<__int32> calculateLimits(const Histogram & hist, const unsigned & delSize)
{
    unsigned lower = std::max(delSize - hist.lowerQuantileDist, delSize / 2);
    return Pair<__int32>(lower, delSize + hist.upperQuantileDist);
}
// Overload for multiple read groups 
inline void calculateLimits(String<Pair<__int32> > & limits,
                            const String<Histogram> & hists,
                            const unsigned & delSize)
{
    resize(limits, length(hists));
    for (unsigned i = 0; i < length(hists); ++i)
    {
        limits[i] = calculateLimits(hists[i], delSize);
    }
}
// ======================================================================================
// Function updateSupportCandidate()
// ======================================================================================
// Check the conditions and update the start- or end position of the candidate if required.
inline void updateSupportCandidate(SupportCandidate & candidate,
                                   const ChromosomeProfile & chromosomeProfile,
                                   const TRGs & rgs,
                                   const String<Histogram> & hists,
                                   const Call & call)
{
    if (checkLRPass(call))
    {
        unsigned tmpDelSize;
        if (candidate.fresh)
        {
            tmpDelSize = call.deletionLength;
        }
        else
        {
            if (candidate.bestLikelihoodRatio < call.likelihoodRatio)
            {
                candidate.frequency = call.frequency;
                candidate.delSize = call.deletionLength;
                candidate.bestLikelihoodRatio = call.likelihoodRatio;
            }
            tmpDelSize = candidate.delSize;
        }
        String<Pair<__int32> > limits;
        calculateLimits(limits, hists, tmpDelSize);
        if (chromosomeProfile.checkActiveReads(rgs, limits))
        {
            if (candidate.start == maxValue<__uint32>())
            {
                candidate.start = call.position;
                candidate.frequency = call.frequency;
                candidate.delSize = call.deletionLength;
                candidate.bestLikelihoodRatio = call.likelihoodRatio;
                candidate.fresh = false;                // Protects the start and frequency values.
            }
            else if (stretchTooLong(chromosomeProfile.currentPos, candidate.start, tmpDelSize))
            {
                candidate.end = chromosomeProfile.currentPos;
                candidate.fresh = true;
            }
        }
        else if (candidate.start != maxValue<__uint32>())
        {
            candidate.end = chromosomeProfile.currentPos;
            candidate.fresh = true;
        }
    }
    else
    {
        if (candidate.start == maxValue<__uint32>() || candidate.fresh)
        {
            return;
        }
        else
        {
            if (stretchTooLong(chromosomeProfile.currentPos, candidate.start, candidate.delSize))
            {
                candidate.end = chromosomeProfile.currentPos;
                candidate.fresh = true;
            }
            else
            {
                String<Pair<__int32> > limits;
                calculateLimits(limits, hists, candidate.delSize);
                if (!chromosomeProfile.checkActiveReads(rgs, limits))
                {
                    candidate.end = chromosomeProfile.currentPos;
                    candidate.fresh = true;
                }
            }
        }
    }
}
// ======================================================================================
// Function getDelWinNum()
// ======================================================================================
// Return the max. number of windows spanned by the deletion.
inline unsigned getDelWinNum(unsigned delSize, unsigned windowSize)
{
    return static_cast<unsigned>(std::ceil(static_cast<double>(delSize) / static_cast<double>(windowSize)));
}
// ======================================================================================
// Function getDelEndWin()
// ======================================================================================
// Return the latest possible end window of the deletion.
inline unsigned getDelEndWin(unsigned delStartWin, unsigned delSize, unsigned windowSize)
{
    return delStartWin + getDelWinNum(delSize, windowSize) - 1;
}
// ======================================================================================
// Function getMinInsertSizeToDel()
// ======================================================================================
// Return the min. insert-size needed to reach and span the deletion.
inline unsigned getMinInsertSizeToDel(const SimulationParameters & params, unsigned pos)
{
    int res = std::ceil(1.5 * params.readSize + params.windowSize * (params.delStartWin - pos));
    return res < 0 ? 0:res;
}
// ======================================================================================
// Function getProbability()
// ======================================================================================
// Return the probability of gaining a read-pair with big enough instert size to reach and span the deletion.
inline double getProbability(const SimulationParameters & params, unsigned pos)
{
    unsigned minSizeToDel = getMinInsertSizeToDel(params, pos);
    if (minSizeToDel <= params.offset)
        return 1.0;
    else if (minSizeToDel - params.offset > length(params.cumulativeDensity) - 1)
        return 0.0;
    else
        return 1 - params.cumulativeDensity[minSizeToDel - params.offset];
}
// ======================================================================================
// Function simulateActiveReads()
// ======================================================================================
// Fill expectedReads with the number of expected active read pairs
inline void simulateActiveReads(String<double>& expectedReads, const SimulationParameters & params)
{
    double avgNewReads = params.avgNewActiveReadsPerWin * params.alleleFrequency;
    resize(expectedReads, 2 * params.flankingWindows + getDelWinNum(params.deletionSize, params.windowSize), 0);
    unsigned deletionEndWin = getDelEndWin(params.delStartWin, params.deletionSize, params.windowSize);
    unsigned i = 0;
    for (; i < params.delStartWin; ++i)
    {
        expectedReads[i] = avgNewReads * getProbability(params, i);
    }
    expectedReads[params.delStartWin] = avgNewReads * 3; // TODO: Fine-tune factor. Depends on w and read-length.
    i = deletionEndWin - 1;
    for (; i < deletionEndWin + params. flankingWindows; ++i)
    {
        unsigned winSinceEnd = i - deletionEndWin;
        expectedReads[i] = -1 * expectedReads[params.delStartWin - winSinceEnd];
    }
    std::partial_sum(begin(expectedReads), end(expectedReads), begin(expectedReads));
    // Fix rounding errors at end of string.
    Iterator<String<double> >::Type it = end(expectedReads);
    while (it > begin(expectedReads))
    {
        --it;
        if (*it < 0.0000000001)
            *it = 0;
        else
            return;
    }
}
// Adds the new simulated active reads to the old set, increasing the count. //TODO: Can be done more efficiently.
inline void simulateActiveReads(String<double>& expectedReads, const SimulationParameters & params, Add add)
{
    (void)add;                                                              // Just touch the tag to prevent warnings.
    String<double> newExpectedReads;
    resize(newExpectedReads, length(expectedReads), 0);
    simulateActiveReads(newExpectedReads, params);
    for(unsigned i = 0; i < length(expectedReads); ++i)
        expectedReads[i] += newExpectedReads[i];
}
// Only simulate up to the delStart window.
inline void simulateActiveReads(String<double>& expectedReads, const SimulationParameters & params, FirstHalf fH)
{
    (void)fH;
    double avgNewReads = params.avgNewActiveReadsPerWin * params.alleleFrequency;
    resize(expectedReads, params.flankingWindows + 1, 0);
    for (unsigned i = 0; i < params.flankingWindows; ++i)
    {
        expectedReads[i] = avgNewReads * getProbability(params, i);
    }
    expectedReads[params.flankingWindows] = avgNewReads * 4; // TODO: Fine-tune factor. Depends on w and read-length.
    //std::partial_sum(begin(expectedReads), end(expectedReads), begin(expectedReads));
    // Fix rounding errors at end of string.
    Iterator<String<double> >::Type it = end(expectedReads);
    while (it > begin(expectedReads))
    {
        --it;
        if (*it < 0.0000000001)
            *it = 0;
        else
            return;
    }
}
// Adds the new simulated active reads to the old set, increasing the count. //TODO: Can be done more efficiently.
// Only performs the simulation up to the assumed start position of the deletion.
inline void simulateActiveReads(String<double>& expectedReads,
                                const SimulationParameters & params,
                                Add add,
                                FirstHalf fH)
{
    (void)add;                                                              // Just touch the tag to prevent warnings.
    String<double> newExpectedReads;
    resize(newExpectedReads, length(expectedReads), 0);
    simulateActiveReads(newExpectedReads, params, fH);
    for(unsigned i = 0; i < length(expectedReads); ++i)
        expectedReads[i] += newExpectedReads[i];
}
// ======================================================================================
// Function calculateSupDist()
// ======================================================================================
// Returns the maximum distance of the cumulative values at the same position between aBegin to aEnd and bBegin to bEnd.
// The number of values in between must be identical. See D in Kolmogorov-Smirnov-test.
inline double calculateSupDist(Iterator<const double>::Type aBegin,
                               Iterator<const double>::Type aEnd,
                               Iterator<const double>::Type bBegin,
                               Iterator<const double>::Type bEnd)
{
    (void)bEnd;                 // Touch it, to suppress warning. (var is only used in debug builds)
    double currentMax = 0;
    double aSum = 0;
    double bSum = 0;
    while(aBegin != aEnd)
    {
        SEQAN_ASSERT_NEQ(bBegin, bEnd);
        aSum += *aBegin;
        bSum += *bBegin;
        currentMax = std::max(std::abs(aSum - bSum), currentMax);
        ++aBegin;
        ++bBegin;
    }
    return currentMax;
}
// ======================================================================================
// Function getBestStartingWindows()
// ======================================================================================
inline void getBestStartingWindows(SupportChangeMap & suppMap,
                                   ChromosomeProfile & chromosomeProfile,
                                   PopDelCallParameters & params,
                                   const unsigned & flankingWindows)
{
    Iterator<String<SupportStretch>, Rooted>::Type suppMapIt = begin(suppMap.supportStretches, Rooted());
    __uint32 maxWinNum = 0;
    for (unsigned i = 0; i < length(suppMap.supportStretches); ++i)
    {
        maxWinNum = std::max(maxWinNum, (suppMap.supportStretches[i].numWindows));
    }
    String<double> activeReads;
    resize(activeReads, maxWinNum + flankingWindows, 0);
    while (!atEnd(suppMapIt))
    {
        for (unsigned i = 0; i < length(activeReads); ++i)
            activeReads[i] = 0;
        unsigned lower;
        unsigned correctedFlankingWindows = flankingWindows;
        if (static_cast<int>(chromosomeProfile.globalMinPos) > static_cast<int>(suppMapIt->start) -
                                                               static_cast<int>(flankingWindows * params.windowShift))
        {
            lower = chromosomeProfile.globalMinPos;
            correctedFlankingWindows = (suppMapIt->start - chromosomeProfile.globalMinPos) / params.windowShift;
        }
        else
        {
            lower = suppMapIt->start - flankingWindows * params.windowShift;
        }
        unsigned upper = suppMapIt->start + std::round(static_cast<double>(suppMapIt->delSize) / 3);  // Limit calculations to ~ start + 1/3 deletion size.
        chromosomeProfile.goToPosition(lower, params.windowShift);
        Iterator<String<double>, Standard>::Type arIt =  begin(activeReads, Standard());
        Iterator<String<double>, Standard>::Type arItEnd =  end(activeReads);
        while(chromosomeProfile.currentPos <= upper && arIt != arItEnd)
        {
            // counting active reads for all RGs in this window.
            for (unsigned rg = 0; rg < length(chromosomeProfile.activeReads); ++rg)
            {
                int minLen= suppMapIt->delSize - 1.96 * params.histograms[rg].stddev;
                *arIt += chromosomeProfile.getActiveReadsNum(rg, chromosomeProfile.currentPos, minLen);
            }
            if (!chromosomeProfile.nextWindow(params.windowShift))
                break;
            ++arIt;
        }
        maxWinNum = suppMapIt->numWindows + correctedFlankingWindows;
        // Create and fill the simulation vector
        String<double> expectedReads;
        resize(expectedReads, correctedFlankingWindows + 1, 0);
        for (unsigned rg = 0; rg < length(params.histograms); ++rg)
        {
            SimulationParameters simParams(params.histograms[rg].cumulativeDensity,
                                           params.windowSize,
                                           suppMapIt->delSize,
                                           params.histograms[rg].readLength,
                                           chromosomeProfile.avgNewReadsPerWindow[rg],
                                           suppMapIt->frequency,
                                           params.histograms[rg].offset,
                                           correctedFlankingWindows);
            simulateActiveReads(expectedReads, simParams, Add(), FirstHalf());
        }
        // Now compare with the simulation for each window.
        Iterator<const String<double>, Standard>::Type exIt = begin(expectedReads, Standard());
        Iterator<const String<double>, Standard>::Type exItEnd = end(expectedReads);
        double currentMinDist = maxValue<double>();
        unsigned bestIndex = 0;
        for (unsigned i = 0; i < std::round(static_cast<double>(suppMapIt->numWindows) / 3); ++i)
        {
            Iterator<const String<double>, Standard>::Type obsIt = begin(activeReads, Standard()) + i;
            Iterator<const String<double>, Standard>::Type obsItEnd = obsIt + (correctedFlankingWindows + 1);
            double dist = calculateSupDist(exIt, exItEnd, obsIt, obsItEnd);
            if (currentMinDist > dist)
            {
                currentMinDist = dist;
                bestIndex = i;
            }
        }
        suppMapIt->bestWindow = suppMapIt->start + bestIndex * params.windowSize;
        ++suppMapIt;
    }
}
#endif /*DISTRIBUTION_MODEL_POPDEL_H_*/
