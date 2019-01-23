#ifndef GENOTYPE_DELETION_POPDEL_CALL_H_
#define GENOTYPE_DELETION_POPDEL_CALL_H_

#include <math.h>
#include "parameter_parsing_popdel_call.h"
#include "profile_structure_popdel_call.h"
#include "active_reads_distribution_popdel_call.h"

using namespace seqan;

// -----------------------------------------------------------------------------
// Function upperHalfMedian()
// -----------------------------------------------------------------------------
// Find and and return the rounded median of the upper half of "values" aka the Q3/third quartile.
template<typename TValue>
inline int upperHalfMedian(String<TValue> & values)
{
    std::sort(begin(values), end(values));
    unsigned n = length(values);
    if (n == 0)
        return 0;
    if (n < 4)
        return _round(values[n - 1]);
    double pos = (3.0 * n + 2.0 + (n % 2)) / 4.0 - 1.0;               // Moore and McCabe (2002) (corrected for 0 index)
    unsigned l = pos;                                                 // Integer part of position
    double r = pos - l;                                               // Rational part of position
    return _round((1 - r) * values[l] + r * values[l + 1]);
}
// -----------------------------------------------------------------------------
// Function initialize_deletion_lengths()
// -----------------------------------------------------------------------------
// Initialize the size of deletions as the upper half median of insert sizes. Use only those which are within
// 50 BP proximity and above the thresholds. Mark the samples as lowCoverage in lowCoverageSample if cov < 2.
inline std::set<int> initialize_deletion_lengths(const ChromosomeProfile & chromosomeProfiles,  //TODO: dont return set.
                                                 const TRGs & rgs,
                                                 const String<unsigned> & thresholds,
                                                 String<bool> & lowCoverageSamples)
{
    unsigned sampleCount = length(rgs);
    String<int> deviations;                                     // Index corresponds to sample
    reserve(deviations, sampleCount);
    // Compute upper half median for each sample
    for (unsigned i = 0; i < sampleCount; ++i)                  // For each sample...
    {
        String<int> sampleValues;
        if(chromosomeProfiles.getActiveReadsDeviations(sampleValues, rgs[i], 2u) < 2)
        {
            lowCoverageSamples[i] = true;
            continue;                                                                    // TODO: Parameter min cov.
        }
        appendValue(deviations, upperHalfMedian(sampleValues));     // TODO: Try different initializations for del-size.
    }
    // Keep only averages of medians that are less than 50 bp apart  and those that are above the threshold.
    std::set<int> deletion_lengths;
    if (empty(deviations))
        return deletion_lengths;
    std::sort(begin(deviations), end(deviations));
    int sum = deviations[0];
    int n = 1;
    unsigned threshold = thresholds[0];
    for (unsigned i = 1; i < length(deviations); ++i)
    {
        if (deviations[i-1] + 50 > deviations[i])
        {
            sum += deviations[i];
            ++n;
            threshold = std::min(threshold, thresholds[i]);
        }
        else
        {
            if (sum/n > (int)threshold)
                deletion_lengths.insert(sum/n);
            sum = deviations[i];
            n = 1;
            threshold = thresholds[i];
        }
    }
    if (sum/n > static_cast<int>(threshold))
        deletion_lengths.insert(sum/n);

    return deletion_lengths;
}
// -----------------------------------------------------------------------------
// Function initialize_allele_frequency()
// -----------------------------------------------------------------------------
// Return a first estimate for the allele frequency based on the number of read pairs whose insert size deviation
// lie in a window of "deletion_length" +- 2x stddev.
inline double initialize_allele_frequency(const ChromosomeProfile & chromosomeProfiles,
                                          const TRGs & rgs,
                                          unsigned deletion_length,
                                          const String<Histogram> & histograms)
{
    unsigned globalCount = 0;
    unsigned globalTotal = 0;
    for (unsigned i = 0; i < length(rgs); ++i)                  // For each sample...
    {
        unsigned count = 0;
        unsigned total = 0;
        const String<unsigned>& currentSample = rgs[i];
        for (unsigned j = 0; j < length(currentSample); ++j)        // For each read group of the current sample...
        {
            __uint32 rg = currentSample[j];
            unsigned currentRgActiveReadCount = chromosomeProfiles.activeReads[rg].size();
            if (currentRgActiveReadCount == 0)
                continue;
            total += currentRgActiveReadCount;
            // Define a window of size 4 * stddev around the deletion length.
            int windowBegin = std::max(static_cast<int>(deletion_length) / 2,
                                       _round(deletion_length - 2 * histograms[rg].stddev));
            int windowEnd = deletion_length + 2 * histograms[rg].stddev;
            ChromosomeProfile::TActiveSet::const_iterator where(chromosomeProfiles.activeReads[rg].begin());
            ChromosomeProfile::TActiveSet::const_iterator whereEnd(chromosomeProfiles.activeReads[rg].end());
            while (where != whereEnd)                             // Count number of  of values(=insert size deviations)
            {                                                     // that fall into the window.
                int deviation = chromosomeProfiles.startProfiles[rg].getDeviationAt(*where);
                if (deviation > windowBegin && deviation < windowEnd)
                    ++count;
                ++where;
            }
        }
        globalCount += count;
        globalTotal += total;
    }
    return static_cast<double>(globalCount) / globalTotal;
}
// -----------------------------------------------------------------------------
// Function assignDad()
// -----------------------------------------------------------------------------
inline void assignDad(Dad & dad,
                      const int & deviation,
                      const int & refUpperBoundary,
                      const int & delLowerBoundary,
                      const int & delUpperBoundary)
{
    if (deviation > refUpperBoundary)
    {
        if (deviation < delLowerBoundary)
        {
            ++dad.between;
        }
        else
        {
            if (deviation <= delUpperBoundary)
            {
                ++dad.alt;
            }
            else
            {
                ++dad.right;
            }
        }
    }
    else
    {
        if (deviation < delUpperBoundary)
        {
            ++dad.ref;
        }
        else
        {
            ++dad.both;
        }
    }
}
// -----------------------------------------------------------------------------
// Function compute_data_likelihoods()
// -----------------------------------------------------------------------------
// Compute the log-likelihoods of the three possible genotypes:
// Homozygous no deletion, heterozygous deletion, homozygous deletion.
// Return a triple of these three log likelihoods scaled by the highest of the three log-likelihoods.
inline Triple<long double> compute_data_likelihoods(String<Triple<long double> > & rgWiseDataLikelihoods,
                                                    const ChromosomeProfile & chromosomeProfiles,
                                                    const String<unsigned> & sample,
                                                    unsigned deletion_length,
                                                    const String<int> & referenceShifts,
                                                    const String<Histogram> & hists)
{
    Triple<long double> logLikelihoods(0, 0, 0);       // log likelihoods for Hom_noDel, Het_del, Hom_del in this order.
    for (unsigned i = 0; i < length(sample); ++i)
    {
        __uint32 rg = sample[i];
        int refShift = referenceShifts[rg];
        Triple<long double> & currentRgWiseDataLikelihoods = rgWiseDataLikelihoods[rg];
        currentRgWiseDataLikelihoods = Triple<long double>(0, 0, 0);
        const Histogram & hist = hists[rg];
        ChromosomeProfile::TActiveSet::const_iterator it(chromosomeProfiles.activeReads[rg].begin());
        ChromosomeProfile::TActiveSet::const_iterator itEnd(chromosomeProfiles.activeReads[rg].end());
        while (it != itEnd)
        {
            int currentDeviation = chromosomeProfiles.getSingleDeviation(rg, it);
            long double g0 = log(I(hist, currentDeviation - refShift));
            long double g1 = log(I(hist, currentDeviation - refShift) +
                                 I(hist, currentDeviation - deletion_length)) -
                             log(2.0);
            long double g2 = log(I(hist, currentDeviation - deletion_length));
            currentRgWiseDataLikelihoods.i1 += g0;
            currentRgWiseDataLikelihoods.i2 += g1;
            currentRgWiseDataLikelihoods.i3 += g2;
            logLikelihoods.i1 += g0;
            logLikelihoods.i2 += g1;
            logLikelihoods.i3 += g2;
            ++it;
        }
        long double maxRgDl = std::max(std::max(currentRgWiseDataLikelihoods.i1, currentRgWiseDataLikelihoods.i2),
                                       currentRgWiseDataLikelihoods.i3);
        currentRgWiseDataLikelihoods.i1 = exp(currentRgWiseDataLikelihoods.i1 - maxRgDl);
        currentRgWiseDataLikelihoods.i2 = exp(currentRgWiseDataLikelihoods.i2 - maxRgDl);
        currentRgWiseDataLikelihoods.i3 = exp(currentRgWiseDataLikelihoods.i3 - maxRgDl);
    }
    // Scale the likelihoods with the largest of the three genotypes.
    long double max_gt = std::max(std::max(logLikelihoods.i1, logLikelihoods.i2), logLikelihoods.i3);
    Triple<long double> res = Triple<long double>(exp(logLikelihoods.i1 - max_gt), //Apply exp to get rid of of the log
                                                  exp(logLikelihoods.i2 - max_gt),
                                                  exp(logLikelihoods.i3 - max_gt));
    // If the likelihoods are too small res will contain 0's. This should only happen for extremely high coverage
    if (res.i1 == 0 || res.i2 == 0 || res.i3 == 0)
    {
        res.i1 = 1;                 // For high coverage regions: Assume reference to avoid calls. TODO: Check
        res.i2 = 0.0000000001;
        res.i3 = 0.0000000001;
    }
    if (res.i1 == res.i2 && res.i1 == res.i3)
    {
        res.i1 = 1;                 // If everything fits equally good/bad assume reference.
        res.i2 = 0.0000000001;
        res.i3 = 0.0000000001;
    }
    return res;
}
//Overload, giving information on the genotype likelihoods and LAD, DAD and mappDist.
inline Triple<long double> compute_data_likelihoods(Triple<long double> & gtLogs,
                                                    Triple<unsigned> & lad,
                                                    Dad & dad,
                                                    Pair<unsigned> & firstLast,
                                                    Pair<unsigned> & suppFirstLast,
                                                    const ChromosomeProfile & chromosomeProfiles,
                                                    const String<unsigned> & sample,
                                                    unsigned deletion_length,
                                                    const String<int> & referenceShifts,
                                                    const String<Histogram> & hists)
{
    Triple<long double> logLikelihoods(0, 0, 0);       // log likelihoods for Hom_noDel, Het_del, Hom_del in this order.
    gtLogs = logLikelihoods;
    int delLowerBorder = maxValue<int>();
    int delUpperBorder = 0;
    for (unsigned i = 0; i < length(sample); ++i)
    {
        __uint32 rg = sample[i];
        const Histogram & hist = hists[rg];
        int refShift = referenceShifts[rg];
        delLowerBorder =  deletion_length - hist.lowerQuantileDist;  //TODO: Might be too strict
        delUpperBorder =  deletion_length + hist.upperQuantileDist;
        ChromosomeProfile::TActiveSet::const_iterator it(chromosomeProfiles.activeReads[rg].begin());
        ChromosomeProfile::TActiveSet::const_iterator itEnd(chromosomeProfiles.activeReads[rg].end());
        while (it != itEnd)
        {
            int currentDeviation = chromosomeProfiles.getSingleDeviation(rg, it);
            assignDad(dad, currentDeviation, hist.upperQuantileDist, delLowerBorder, delUpperBorder);
            double refLikelihood = I(hist, currentDeviation - refShift);
            double delLikelihood = I(hist, currentDeviation - deletion_length);
            if (refLikelihood >= 2 * delLikelihood)
                ++lad.i1;       // Supporting the reference
            else if (delLikelihood >= 2 * refLikelihood)
                ++lad.i3;       // Supporting a deletion
            else
                ++lad.i2;       // Unclear support.
            logLikelihoods.i1 += log(refLikelihood);
            gtLogs.i1 += log10(I(hist, currentDeviation - refShift));
            logLikelihoods.i2 += log(I(hist, currentDeviation - refShift) +
                                     I(hist, currentDeviation - deletion_length)) -
                                 log(2.0);
            gtLogs.i2 += log10(I(hist, currentDeviation - refShift) +
                               I(hist, currentDeviation - deletion_length)) -
                         log10(2.0);
            logLikelihoods.i3 += log(delLikelihood);
            gtLogs.i3 += log10(I(hist, currentDeviation - deletion_length));
            ++it;
        }
    }
    firstLast = chromosomeProfiles.getActiveReadsFirstLast(hists, sample);
    chromosomeProfiles.updateSupportFirstLast(suppFirstLast, hists, sample, delLowerBorder, delUpperBorder);
    // Scale the likelihoods with the largest of the three genotypes.
    long double max_gt = std::max(std::max(gtLogs.i1, gtLogs.i2), gtLogs.i3);
    long double max_dl = std::max(std::max(logLikelihoods.i1, logLikelihoods.i2), logLikelihoods.i3);
    gtLogs.i1 -= max_gt;
    gtLogs.i2 -= max_gt;
    gtLogs.i3 -= max_gt;
    if (gtLogs.i1  == gtLogs.i2  && gtLogs.i1 == gtLogs.i3) // TODO: Check if beneficial for nested deletions.
    {
        gtLogs.i1 = 0;
        gtLogs.i2 = -10;
        gtLogs.i3 = -10;
    }
    Triple<long double> res = Triple<long double>(exp(logLikelihoods.i1 - max_dl), //Apply exp to get rid of of the log
                                                  exp(logLikelihoods.i2 - max_dl),
                                                  exp(logLikelihoods.i3 - max_dl));
    // If the likelihoods are too small res will contain 0's. This should only happen for extremely high coverage
    if (res.i1 == 0 || res.i2 == 0 || res.i3 == 0)
    {
        res.i1 = 1;                 // For high coverage regions: Assume reference to avoid calls. TODO: Check
        res.i2 = 0.0000000001;
        res.i3 = 0.0000000001;
    }
    if (res.i1 == res.i2 && res.i1 == res.i3)    // TODO: Check if beneficial for nested deletions.
    {
        res.i1 = 1;                 // If everything fits equally good/bad assume reference.
        res.i2 = 0.0000000001;
        res.i3 = 0.0000000001;
    }
    return res;
}
// -----------------------------------------------------------------------------
// Function compute_gt_likelihoods()
// -----------------------------------------------------------------------------
// Return the computed genotype likelihoods given the allele frequency.
inline Triple<double> compute_gt_likelihoods(double allele_frequency)
{
    SEQAN_ASSERT_GEQ(allele_frequency, 0.0);
    double pseudoFreq = 0.0000000001;
    Triple<double> gt_likelihoods;
    gt_likelihoods.i1 = std::max((1 - allele_frequency) * (1 - allele_frequency), pseudoFreq);
    gt_likelihoods.i2 = std::max(2 * allele_frequency * (1 - allele_frequency), pseudoFreq);
    gt_likelihoods.i3 = std::max(allele_frequency * allele_frequency, pseudoFreq);
    return gt_likelihoods;
}
// -----------------------------------------------------------------------------
// Function update_deletion_length()    //TODO: Export RefShift into own function.
// -----------------------------------------------------------------------------
// Update the estimate for the deletion length considering the data and GT likelihood of each sample and read group.
// Return Pair<updated reference shift, updated deletion length>
inline unsigned update_deletion_length(const ChromosomeProfile & chromosomeProfiles,
                                       const TRGs & rgs,
                                       const String<Triple<long double> > & data_likelihoods,
                                       const String<Triple<long double> > & rgWiseDataLikelihoods,
                                       const Triple<double> & gt_likelihoods,
                                       int deletion_length,
                                       String<int> & referenceShifts,
                                       const String<Histogram> & histograms)
{
    // Compute the sum and weighted sum of probabilities that the read pairs are from a deletion.
    double sumDel = 0;
    double weighted_sumDel = 0;
    Iterator<const String<Triple<long double> > >::Type dlIt = begin(data_likelihoods);
    Iterator<const String<Triple<long double> > >::Type rgDlIt = begin(rgWiseDataLikelihoods);
    Iterator<const String<String<unsigned> > >::Type sIt = begin(rgs);
    Iterator<const String<String<unsigned> > >::Type sItEnd = end(rgs);
    while (sIt != sItEnd)
    {

        double aSum = getValueI1(*dlIt) * gt_likelihoods.i1 +
                      getValueI2(*dlIt) * gt_likelihoods.i2 +
                      getValueI3(*dlIt) * gt_likelihoods.i3;
        double a1 = log(getValueI2(*dlIt)) + log(gt_likelihoods.i2) - log(aSum);
        double a2 = log(getValueI3(*dlIt)) + log(gt_likelihoods.i3) - log(aSum);
        ++dlIt;
        Iterator <const TReadGroupIndices>::Type rgIt = begin(*sIt);
        Iterator <const TReadGroupIndices>::Type rgItEnd = end(*sIt);
        while (rgIt != rgItEnd)
        {
            double sumRef = 0;
            double weighted_sumRef = 0;
            double aSumRg = getValueI1(*rgDlIt) * gt_likelihoods.i1 +
                            getValueI2(*rgDlIt) * gt_likelihoods.i2 +
                            getValueI3(*rgDlIt) * gt_likelihoods.i3;
            long double l = log(getValueI1(*rgDlIt));
            long double g = log(gt_likelihoods.i1);
            long double a = log(aSumRg);
            long double a0Rg = l + g - a;
            long double a1Rg = log(getValueI2(*rgDlIt)) + log(gt_likelihoods.i2) - log(aSumRg);
            ChromosomeProfile::TActiveSet::const_iterator actReadsIt(chromosomeProfiles.activeReads[*rgIt].begin());
            ChromosomeProfile::TActiveSet::const_iterator actReadsItEnd(chromosomeProfiles.activeReads[*rgIt].end());
            while (actReadsIt != actReadsItEnd)
            {
                int currentDeviation = chromosomeProfiles.getSingleDeviation(*rgIt, actReadsIt);
                double del = I(histograms[*rgIt], currentDeviation - deletion_length);
                double no_del = I(histograms[*rgIt], currentDeviation - referenceShifts[*rgIt]);
                double prob_del = exp(a1)   * del / (del + no_del) + exp(a2);
                double prob_ref = exp(a1Rg) * no_del / (del + no_del) + exp(a0Rg);
                sumDel += prob_del;
                sumRef += prob_ref;
                weighted_sumDel += prob_del * (currentDeviation);
                weighted_sumRef += prob_ref * (currentDeviation);
                ++actReadsIt;
            }
            referenceShifts[*rgIt] = weighted_sumRef / sumRef;
            if (referenceShifts[*rgIt] > histograms[*rgIt].stddev ||
                referenceShifts[*rgIt] < -1 * histograms[*rgIt].stddev)
                referenceShifts[*rgIt] = 0;     // Limits the refShift to 0 if it is too big.
            ++rgIt;
        }
        ++sIt;
    }
    double len = (weighted_sumDel / sumDel);
    if (len < 0)
        return 0;
    else
        return round(len);
}
// -----------------------------------------------------------------------------
// Function update_allele_frequency()
// -----------------------------------------------------------------------------
// Return the updated allele frequency using the data- and genotype-likelihoods.
inline double update_allele_frequency(const String<Triple<long double> > & data_likelihoods,
                                      const Triple<double> & gt_likelihoods)
{
    typedef Iterator<const String<Triple<long double> > >::Type TIter;
    double sum = 0;
    TIter it = begin(data_likelihoods);
    TIter itEnd = end(data_likelihoods);
    while (it != itEnd)
    {
        double p0 = (*it).i1 * gt_likelihoods.i1;
        double p1 = (*it).i2 * gt_likelihoods.i2;
        double p2 = (*it).i3 * gt_likelihoods.i3;
        double pAll = p0 + p1 + p2;
        sum += (p1 + 2 * p2) / pAll;
        ++it;
    }
    SEQAN_ASSERT_GEQ(sum / 2.0 / length(data_likelihoods), 0);
    return sum / 2.0 / length(data_likelihoods);
}
// -----------------------------------------------------------------------------
// Function deletion_likelihood_ratio()
// -----------------------------------------------------------------------------
// Calculates and returns the loglikelihood ratio log(L(del)/L(noDel)).
inline double deletion_likelihood_ratio(const String<Triple<long double> > & data_likelihoods,
                                        const Triple<double> & gt_likelihoods)
{
    double del = 0;
    double no_del = 0;
    for (unsigned i = 0; i < length(data_likelihoods); ++i)
    {
        double p0 = data_likelihoods[i].i1 * gt_likelihoods.i1;
        double p1 = data_likelihoods[i].i2 * gt_likelihoods.i2;
        double p2 = data_likelihoods[i].i3 * gt_likelihoods.i3;
        double pAll = p0 + p1 + p2;
        double a0 = p0 / pAll;
        double a1 = p1 / pAll;
        double a2 = p2 / pAll;
        del += log(a0 * data_likelihoods[i].i1 + a1 * data_likelihoods[i].i2 + a2 * data_likelihoods[i].i3);
        no_del += log(data_likelihoods[i].i1);
    }
    return del - no_del;
}
// -----------------------------------------------------------------------------
// Function genotype_deletion_window()
// -----------------------------------------------------------------------------
// Perform the estimation of the deletion length, allele frequency and likelihood calculations for all samples
// in the window and perform the genotype estimation based on these values.
// Return true if at least one deletion is called (not neccessarily passed).
inline bool genotype_deletion_window(String<Call> & calls,
                                     const ChromosomeProfile & chromosomeProfiles,
                                     const TRGs & rgs,
                                     PopDelCallParameters & params)
{
    // Initialize the deletion length and check if the samples have sufficient coverage.
    String<bool> lowCoverageSamples;
    resize(lowCoverageSamples, length(rgs), false);
    std::set<int> deletion_lengths = initialize_deletion_lengths(chromosomeProfiles,
                                                                 rgs,
                                                                 params.minInitDelLengths,
                                                                 lowCoverageSamples);
    if (deletion_lengths.empty())
        return false;

    bool ret = false;
    String<int> referenceShifts;
    resize(referenceShifts, length(params.histograms));
    for (std::set<int>::iterator it = deletion_lengths.begin(); it != deletion_lengths.end(); ++it)
    {
        for (auto refIt = begin(referenceShifts); refIt != end(referenceShifts); ++refIt)
            *refIt = 0;

        // Keep track of the deletion length estimates in order to notice convergence.
        std::map<int, double> visited;
        // Initialize the allele frequency.
        double freq = initialize_allele_frequency(chromosomeProfiles, rgs, *it, params.histograms);

        if (freq == 0)
            continue;

        // Initialize likelihoods of data given genotype.
        String<Triple<long double> > rgWiseDataLikelihoods;
        resize(rgWiseDataLikelihoods, length(params.histograms));
        String<Triple<long double> > data_likelihoods;
        resize(data_likelihoods, length(rgs));
        for (unsigned i = 0; i < length(rgs); ++i)
            data_likelihoods[i] = compute_data_likelihoods(rgWiseDataLikelihoods,
                                                           chromosomeProfiles,
                                                           rgs[i],
                                                           *it,
                                                           referenceShifts,
                                                           params.histograms);
        // Compute likelihood of genotypes given the allele frequency
        Triple<double> gtLikelihoods = compute_gt_likelihoods(freq);
        // Update the deletion length using the frequency estimate.
        unsigned len = *it;

        double prevFreq = freq;
        unsigned prevLen = len;
        Triple<double> prevGtLikelihoods = gtLikelihoods;
        String<int> prevShifts = referenceShifts;

        __uint32 iterations = 0;
        // Alternate between updating the allele freq and the deletion length until convergence of the deletion length.
        while (len >= params.minLen && iterations < params.iterations)
        {
            ++iterations;
            visited[len] = freq;
            prevLen = len;
            prevFreq = freq;
            prevGtLikelihoods = gtLikelihoods;

            len = update_deletion_length(chromosomeProfiles,
                                         rgs,
                                         data_likelihoods,
                                         rgWiseDataLikelihoods,
                                         gtLikelihoods,
                                         len,
                                         referenceShifts,
                                         params.histograms);

            for (unsigned i = 0; i < length(rgs); ++i)
                data_likelihoods[i] = compute_data_likelihoods(rgWiseDataLikelihoods,
                                                               chromosomeProfiles,
                                                               rgs[i],
                                                               len,
                                                               referenceShifts,
                                                               params.histograms);

            freq = update_allele_frequency(data_likelihoods, gtLikelihoods);
            if (freq == 0)
                break;
            gtLikelihoods = compute_gt_likelihoods(freq);

            // Check if this pair of deletion length and allele frequency has already been observed. Break if TRUE.
            std::map<int,double>::iterator v = visited.find(len);
            if (v != visited.end())
            {
                if (fabs(v->second - freq) <= 0.0001)
                {
                    // Calculate likelihood ratio for the current values.
                    double logLR = deletion_likelihood_ratio(data_likelihoods, gtLikelihoods);

                    // Re-calculate data likelihood for previous values and calculate the likelihood ratio.
                    for (unsigned i = 0; i < length(rgs); ++i)
                        data_likelihoods[i] = compute_data_likelihoods(rgWiseDataLikelihoods,
                                                                       chromosomeProfiles,
                                                                       rgs[i],
                                                                       prevLen,
                                                                       prevShifts,
                                                                       params.histograms);

                    double prevLogLR = deletion_likelihood_ratio(data_likelihoods, prevGtLikelihoods);
                    // Keep better of both estimates.
                    if (prevLogLR > logLR)
                    {
                        len = prevLen;
                        freq = prevFreq;
                        referenceShifts = prevShifts;
                    }
                    break;
                }
            }
        }
        if (freq < 0.0000000001 || len < params.minLen)
        {
            continue;
        }
        String<Triple<long double> > gtLogs;
        resize (gtLogs, length(rgs));
        String<Triple<unsigned> > lads; // Likelihhod based counts of read pairs for the classes.
        resize (lads, length(rgs), Triple<unsigned>(0, 0, 0));
        String<Dad> dads;                  // Distribution based counts of RP for the classes.
        resize (dads, length(rgs));
        String<Pair<unsigned> > firstLasts; // Lowest firstWin and highest lastWin for each sample.
        resize (firstLasts, length(rgs));
        Pair<unsigned> suppFirstLast(0, maxValue<unsigned>());
        for (unsigned i = 0; i < length(rgs); ++i)
            data_likelihoods[i] = compute_data_likelihoods(gtLogs[i],
                                                           lads[i],
                                                           dads[i],
                                                           firstLasts[i],
                                                           suppFirstLast,
                                                           chromosomeProfiles,
                                                           rgs[i],
                                                           len,
                                                           referenceShifts,
                                                           params.histograms);

        if (suppFirstLast.i1 == 0)
            continue;                   // We don't want calls that don't match their supporting reads.

        double logLR = deletion_likelihood_ratio(data_likelihoods, gtLikelihoods);
        if (logLR >= params.minimumLikelihoodRatio)
        {
            if (params.windowWiseOutput)
            {
                appendValue(calls, Call(*it,
                                        iterations,
                                        len,
                                        logLR,
                                        freq,
                                        chromosomeProfiles.currentPos - 1,
                                        0));
            }
            else
            {
                appendValue(calls, Call(*it,
                                        iterations,
                                        len,
                                        logLR,
                                        freq,
                                        suppFirstLast.i1,
                                        suppFirstLast.i2));
            }
            Call & currentCall = calls[length(calls) - 1];
            set(currentCall.lads, lads);
            set(currentCall.dads, dads);
            set(currentCall.firstLast, firstLasts);
            calculatePhredGL(currentCall, gtLogs);
            ret = true;

            resetFilters(currentCall);
            setNoDataSamples(currentCall, lowCoverageSamples);
            if (!checkSampleNumber(lowCoverageSamples, params.minSampleFraction))
                setSampleFilter(currentCall);
        }
    }
    return ret;
}
#endif /* GENOTYPE_DELETION_POPDEL_CALL_H_ */
