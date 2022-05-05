#ifndef DETECT_JUNCTION_POPDEL_CALL_H_
#define DETECT_JUNCTION_POPDEL_CALL_H_

#include "parameter_parsing_popdel_call.h"
#include "profile_structure_popdel_call.h"

using namespace seqan;

typedef std::set<Iterator<const String<TranslocationWindowEntry> >::Type> TTranslocCluster;
typedef std::map<uint32_t, Pair<TTranslocCluster, bool> >TTranslocClusterMap;  // bool at i2 indicates passed filters
struct TranslocationClusterMap
{
    TTranslocClusterMap clusterMap;

    // Insert the iterator it in the map for the given targetChrom.
    // If no set for targetChrom exists, create one.
    inline void insert(uint32_t targetChrom, Iterator<const String<TranslocationWindowEntry> >::Type it)
    {
        auto clusterIt = clusterMap.find(targetChrom);
        if (clusterIt != clusterMap.end())
        {
            clusterIt->second.i1.insert(it);
        }
        else
        {
            clusterMap[targetChrom].i2 = false;
            clusterMap[targetChrom].i1.insert(it);
        }
    }
    // Return the number of cluster in the map
    inline unsigned size() const
    {
        return clusterMap.size();
    }
    // Return an iterator pointing to the first element of TranslocationClusterMap::clusterMap
    inline TTranslocClusterMap::iterator begin()
    {
        return clusterMap.begin();
    }
    // Return an const_iterator pointing to the first element of TranslocationClusterMap::clusterMap
    inline TTranslocClusterMap::const_iterator begin() const
    {
        return clusterMap.begin();
    }
    // Return an iterator pointing one past the last element of TranslocationClusterMap::clusterMap
    inline TTranslocClusterMap::const_iterator end() const
    {
        return clusterMap.end();
    }
    // Return the number of clusters that have the filter set to pass
    inline unsigned getClusterPasses() const
    {
        unsigned c = 0;
        for (auto it = clusterMap.begin(); it != clusterMap.end(); ++it)
        {
            if (it->second.i2)
                ++c;
        }
        return c;
    }
    // Set all passes to false
    inline void markAsBadSample()
    {
        for (auto it = clusterMap.begin(); it != clusterMap.end(); ++it)
            it->second.i2 = false;
    }
};
// -----------------------------------------------------------------------------
// Function getClusterSize()
// -----------------------------------------------------------------------------
// Return the size of the cluster of the iterator it
inline unsigned getClusterSize(const TTranslocClusterMap::const_iterator it)
{
    return it->second.i1.size();
}
// Return the size of the cluster of the iterator it for the given orientation
inline unsigned getClusterSize(const TTranslocClusterMap::const_iterator it, const Orientation o)
{
    unsigned c = 0;
    for (auto entry = it->second.i1.begin(); entry != it->second.i1.end(); ++entry)
    {
        if ((*entry)->orientation == o)
            ++c;
    }
    return c;
}
// Comparison function for two TTranslocCluster iterators, comparing their mate positions.
inline bool lowerMatePos(const TTranslocCluster::iterator l, const TTranslocCluster::iterator r)
{
    return (*l)->second.pos < (*r)->second.pos;
}
// -----------------------------------------------------------------------------
// Function getClusterSizeByOrientation()
// -----------------------------------------------------------------------------
// Return a pair of counts <FR-Cluster-size, RF-Cluster-size>
// inline Pair<unsigned> getClusterSizeByOrientation(const TTranslocClusterMap::iterator it)
// {
//     unsigned fr = 0;
//     unsigned rf = 0;
//     for (auto entry = it->second.i1.begin(); entry != it->second.i1.end(); ++entry)
//     {
//         if ((*entry)->orientation == Orientation::FR)
//             ++fr;
//         else if ((*entry)->orientation == Orientation::RF)
//             ++rf;
//     }
//     return Pair<unsigned>(fr, rf);
// }
// -----------------------------------------------------------------------------
// Function setClusterFilter()
// -----------------------------------------------------------------------------
// Set the filter of the cluster it is pointing to the value in pass
inline void setClusterFilter(const TTranslocClusterMap::iterator it, const bool pass)
{
    it->second.i2 = pass;
}
// -----------------------------------------------------------------------------
// Function checkSingleOrientationCount()
// -----------------------------------------------------------------------------
// Compare the count of oriented read pairs against the thresholds.
// Return true of both thresholds are passed, false otherwise.
inline bool checkSingleOrientationCount(const unsigned & count,
                                        const unsigned & totalActiveReads,
                                        const unsigned & absThreshold = 2,
                                        const double & relativeThreshold = 0.1)
{
    SEQAN_ASSERT_GEQ(relativeThreshold, 0.0);
    SEQAN_ASSERT_LEQ(relativeThreshold, 1.0);
    SEQAN_ASSERT_LEQ(count, totalActiveReads);
    if (totalActiveReads == 0u)
        return false;
    else
        return (count >= absThreshold) && (static_cast<double>(count) / totalActiveReads >= relativeThreshold);
}
// -----------------------------------------------------------------------------
// Function checkOrientationCounts()
// -----------------------------------------------------------------------------
// Return a string of bools, each indicating if the respective counts of RF/FF/RR read pairs of THIS ONE SAMPLE
// passed the thresholds.
// Note: Does not consider FR read pairs, as they are uninteresing for our break point detection.
inline String<bool> checkOrientationCounts(const OrientationCounts & orientationCounts,
                                           const unsigned absThreshold = 2,
                                           const double relativeThreshold = 0.1)
{
    String<bool> passedOrientations;
    resize(passedOrientations, 3, Exact());
    const unsigned totalActiveReads = sum(orientationCounts);
    passedOrientations[0] = checkSingleOrientationCount(orientationCounts.rfCount,
                                                        totalActiveReads,
                                                        absThreshold,
                                                        relativeThreshold);
    passedOrientations[1] = checkSingleOrientationCount(orientationCounts.ffCount,
                                                        totalActiveReads,
                                                        absThreshold,
                                                        relativeThreshold);
    passedOrientations[2] = checkSingleOrientationCount(orientationCounts.rrCount,
                                                        totalActiveReads,
                                                        absThreshold,
                                                        relativeThreshold);
    return passedOrientations;
}
// -----------------------------------------------------------------------------
// Function clusterTranslocationRecords()
// -----------------------------------------------------------------------------
// Cluster the translocation window entries of the given read groups by chromsome
inline void clusterTranslocationRecords(TranslocationClusterMap & clusteredRecs,
                                        const TranslocationProfile & translocProfiles,
                                        const String<uint32_t> & rgs)
{
    for (auto rgIt = begin(rgs); rgIt != end(rgs); ++rgIt)    // For each read group of sample s.
    {   // Iterate over all records in the window
        for (auto it = translocProfiles.activeReads[*rgIt].i1; it !=  translocProfiles.activeReads[*rgIt].i2; ++it)
        {
            clusteredRecs.insert(it->second.refID, it);
        }
    }
}
// -----------------------------------------------------------------------------
// Function getNumPositionClusters()
// -----------------------------------------------------------------------------
// Take a translocation cluster and return how many sub-clusters are created by the positions of the mates.
// NOTE: Only considers records with the given orientation
inline unsigned getNumPositionClusters(TTranslocClusterMap::const_iterator clusterIt,
                                       const double stddev,
                                       const Orientation o,
                                       const double factor = 3)
{
    auto cluster = clusterIt->second.i1;    // The cluster clusterIt is pointing to
    unsigned clusterSize = cluster.size();
    if (clusterSize < 2)
        return clusterSize;

    unsigned margin = std::ceil(factor * stddev);
    String<unsigned> records;
    reserve(records, clusterSize, Exact());
    for (auto it = cluster.begin(); it != cluster.end(); ++it)
        if ((*it)->orientation == o)
            appendValue(records, (*it)->second.pos);

    if (empty(records))
        return 0;

    std::sort(begin(records), end(records));
    unsigned anchorPos = records[0];
    unsigned c = 1;
    for (auto it = begin(records) + 1; it != end(records); ++it)
    {
        SEQAN_ASSERT_LEQ(anchorPos, *it);
        if (anchorPos + margin < *it)
        {
            anchorPos = *it;
            ++c;
        }
    }
    return c;
}
// -----------------------------------------------------------------------------
// Function getMaxPosClusterSize()
// -----------------------------------------------------------------------------
// Return the size of the biggest positional cluster.
// Store the positions (on the remote chrom!) of the biggest cluster in the bestRemotePos String.
// NOTE: Only considers records with the given orientation
inline unsigned getMaxPosCluster(String<unsigned> & bestLocalPos,
                                 String<unsigned> & bestRemotePos,
                                 TTranslocClusterMap::const_iterator clusterIt,
                                 double stddev,
                                 const Orientation o,
                                 double factor = 3)
{
    auto cluster = clusterIt->second.i1;    // The cluster clusterIt is pointing to
    unsigned clusterSize = cluster.size();
    if (clusterSize < 2)
        return clusterSize;

    unsigned margin = std::ceil(factor * stddev);
    String<TTranslocCluster::iterator> records;
    reserve(records, clusterSize, Exact());
    for (auto it = cluster.begin(); it != cluster.end(); ++it)
        if ((*it)->orientation == o)
            appendValue(records, it); // (*it)->second.pos

    if (empty(records))
        return 0;

    std::sort(begin(records), end(records), lowerMatePos);
    unsigned anchorPos = (*records[0])->second.pos;
    std::multiset<Pair<unsigned> > clusterMax;         // i1: pos on remote chrom. i2. pos on current chrom
    std::multiset<Pair<unsigned> > clusterCurrent;     // i1: pos on remote chrom. i2. pos on current chrom
    clusterCurrent.insert(Pair<unsigned>(anchorPos, (*records[0])->first.pos));
    for (auto it = begin(records) + 1; it != end(records); ++it)
    {
        unsigned localPos  = (**it)->first.pos;
        unsigned remotePos = (**it)->second.pos;
        SEQAN_ASSERT_LEQ(anchorPos, remotePos);
        if (anchorPos + margin < remotePos)
        {
            if (clusterCurrent.size() > clusterMax.size())
            {
                clusterMax = clusterCurrent;
                clusterCurrent.clear();
            }
            anchorPos = remotePos;
            clusterCurrent.insert(Pair<unsigned>(remotePos, localPos));
        }
        else
        {
            clusterCurrent.insert(Pair<unsigned>(remotePos, localPos));
        }
    }
    // Copy contents of biggest cluster to the bestRemotePos String
    if (clusterCurrent.size() > clusterMax.size())
    {
        clusterMax.clear();
        reserve(bestRemotePos, clusterCurrent.size(), Exact());
        reserve(bestLocalPos, clusterCurrent.size(), Exact());
        for (auto it = clusterCurrent.begin(); it != clusterCurrent.end(); ++it)
        {
            appendValue(bestRemotePos, it->i1);
            appendValue(bestLocalPos, it->i2);
        }
    }
    else
    {
        clusterCurrent.clear();
        reserve(bestRemotePos, clusterMax.size(), Exact());
        reserve(bestLocalPos, clusterMax.size(), Exact());
        for (auto it = clusterMax.begin(); it != clusterMax.end(); ++it)
        {
            appendValue(bestRemotePos, it->i1);
            appendValue(bestLocalPos, it->i2);
        }
    }
    // bestRemotePos is already sorted because we used a multiset for storing the positions before.
    std::sort(begin(bestLocalPos), end(bestLocalPos));
    return length(bestRemotePos);
}
// -----------------------------------------------------------------------------
// Function getPassingTargetChroms()
// -----------------------------------------------------------------------------
// Take the translocationMaps and create a set of all target chromosomes from clusters marked as passing.
inline void getPassingTargetChroms(std::set<unsigned> & targetChroms,
                                   const String<TranslocationClusterMap> & clusterMaps,
                                   const String<bool> & sampleWiseTranslocationPasses)
{
    targetChroms.clear();
    for (unsigned s = 0; s < length(clusterMaps); ++s)
    {
        if (!sampleWiseTranslocationPasses[s])  // Skip this sample if it does not have any passes at all.
            continue;

        for (auto it = clusterMaps[s].begin(); it != clusterMaps[s].end(); ++it)
            if (it->second.i2)  // Check if the cluster is marked as passing
                targetChroms.insert(it->first);   // insert the chromosome number into the set.
    }
}
// -----------------------------------------------------------------------------
// Function getLocalAndRemotePos()
// -----------------------------------------------------------------------------
// Take the string of <localPos, remotePos> Pairs (one per sampe) between start, end for both FR and RF and return the
// local and remote breakpoint positions. i1: local, i2. remote
// The final positions are calcualted as the weighted means of the fr and rf clusters
inline Pair<unsigned> getLocalAndRemotePos(const Iterator<String<Pair<unsigned> > >::Type frStart,
                                           const Iterator<String<Pair<unsigned> > >::Type frEnd,
                                           const Iterator<String<Pair<unsigned> > >::Type rfStart,
                                           const Iterator<String<Pair<unsigned> > >::Type rfEnd)
{
    SEQAN_ASSERT_EQ(difference(frStart, frEnd), difference(frStart, frEnd));
    String<unsigned> frLoc;
    String<unsigned> frRem;
    String<unsigned> rfLoc;
    String<unsigned> rfRem;
    reserve(frLoc, difference(frStart, frEnd), Exact());
    reserve(frRem, difference(frStart, frEnd), Exact());
    reserve(rfLoc, difference(rfStart, rfEnd), Exact());
    reserve(rfRem, difference(rfStart, rfEnd), Exact());

    auto frIt = frStart;
    auto rfIt = rfStart;
    while (frIt != frEnd)
    {
        if (frIt->i1 != 0)  // fr[s].i2 should automatically be != 0 if fr[s].i1 != 0
        {
            SEQAN_ASSERT(frIt->i2 != 0u);
            appendValue(frLoc, frIt->i1);
            appendValue(frRem, frIt->i2);
        }
        ++frIt;
    }
    while (rfIt != rfEnd)
    {
        if (rfIt->i1 != 0)
        {
            SEQAN_ASSERT(rfIt->i2 != 0u);
            appendValue(rfLoc, rfIt->i1);
            appendValue(rfRem, rfIt->i2);
        }
        ++rfIt;
    }
    unsigned medianFRloc = 0;
    unsigned medianFRrem = 0;
    unsigned medianRFloc = 0;
    unsigned medianRFrem = 0;
    unsigned frCount = length(frLoc);
    unsigned rfCount = length(rfLoc);
    if (frCount != 0)
    {
        medianFRloc = getPercentile(frLoc, 0.5);
        medianFRrem = getPercentile(frRem, 0.5);
    }
     if (rfCount != 0)
    {
        medianRFloc = getPercentile(rfLoc, 0.5);
        medianRFrem = getPercentile(rfRem, 0.5);
    }
    return Pair<unsigned>((medianFRloc * frCount + medianRFloc * rfCount) / (frCount + rfCount),
                          (medianFRrem * frCount + medianRFrem * rfCount) / (frCount + rfCount));
}
// Wrapper for working on the whole fr and rf strings
inline Pair<unsigned> getLocalAndRemotePos(String<Pair<unsigned> > & fr, String<Pair<unsigned> > & rf)
{
    return getLocalAndRemotePos(begin(fr), end(fr), begin(rf), end(rf));
}
// -----------------------------------------------------------------------------
// Function assignTranslocOrientation()
// -----------------------------------------------------------------------------
// Assign a oientation to the JunctionCall object based on the most mumerous orientations in the clusters.
inline void assignTranslocOrientation(JunctionCall & call, const Tuple<unsigned, 4>  orientationCounter)
{
    if  (orientationCounter.i[0] >= orientationCounter.i[1])
    {
        if (orientationCounter.i[0] >= orientationCounter.i[2])
        {
            if (orientationCounter.i[0] >= orientationCounter.i[3])
                call.orientation = Orientation::FR;
            else
                call.orientation = Orientation::RR;
        }
        else
        {
            if (orientationCounter.i[2] >= orientationCounter.i[3])
                call.orientation = Orientation::FF;
            else
                call.orientation = Orientation::RR;
        }
    }
    else if (orientationCounter.i[1] >= orientationCounter.i[2])
    {
        if (orientationCounter.i[1] >= orientationCounter.i[3])
            call.orientation = Orientation::RF;
        else
            call.orientation = Orientation::RR;
    }
    else if (orientationCounter.i[2] >= orientationCounter.i[3])
    {
        call.orientation = Orientation::FF;
    }
    else
    {
        call.orientation = Orientation::RR;
    }
}
// Overload for counting the orientations manually
inline void assignTranslocOrientation(JunctionCall & call)
{
    Tuple<unsigned, 4> orientationCounter;
    for (unsigned i = 0; i < 4; ++i)
        orientationCounter.i[i] = 0;

    for  (unsigned s = 0; s < length(call.perSampleTranslocSupport); ++s)
        for (unsigned j = 0; j < 4; ++j)
            orientationCounter.i[j] += call.perSampleTranslocSupport[s].i[j];

    assignTranslocOrientation(call, orientationCounter);
}
// -----------------------------------------------------------------------------
// Function processTranslocationCluster()
// -----------------------------------------------------------------------------
inline bool processTranslocationClusters(String<JunctionCall> & calls,
                                         const String<TranslocationClusterMap> & clusterMaps,
                                         const String<OrientationCounts> & sampleWiseOrientations,
                                         const TRGs & rgs,
                                         const String<bool> & sampleWiseTranslocationPasses,
                                         const bool anyTranslocationPass,
                                         const unsigned currentPos,
                                         const PopDelCallParameters & params)
{
    // TODO: Continue here and check if adding FF and RR has broken anythingmmmake 
    if (!anyTranslocationPass)
        return false;

    const unsigned sampleNum = length(rgs);
    std::set<unsigned> targetChroms;
    getPassingTargetChroms(targetChroms, clusterMaps, sampleWiseTranslocationPasses);

    // Create one JunctionCall object for each target chromosome and gather the info from all samples if they are not
    // marked as bad
    for (auto chromIt = targetChroms.begin(); chromIt != targetChroms.end(); ++chromIt)
    {
        JunctionCall call(currentPos);
        call.svtype = SVType::TRL;
        call.mateChromosome = *chromIt;
        Tuple<double, 5> tupleInit;
        for (unsigned i = 0; i < 5; ++i)
            tupleInit.i[i] = 0.0;
        resize(call.perSampleTranslocSupport, sampleNum, tupleInit, Exact());
        resize(call.perSampleFRPositions, sampleNum, Pair<unsigned>(0, 0), Exact());
        resize(call.perSampleRFPositions, sampleNum, Pair<unsigned>(0, 0), Exact());
        resize(call.perSampleFFPositions, sampleNum, Pair<unsigned>(0, 0), Exact());
        resize(call.perSampleRRPositions, sampleNum, Pair<unsigned>(0, 0), Exact());
        Tuple<unsigned, 4> orientationCounter;
        for (unsigned i = 0; i < 4; ++i)
            orientationCounter.i[i] = 0;

        unsigned totalFrRfEvidence = 0;
        unsigned totalFfRrEvidence = 0;
        for (unsigned s = 0; s < sampleNum; ++s)
        {
            if (sampleWiseTranslocationPasses[s])   // Pass for at least one translocation / is not bad
            {
                double stddev = params.histograms[rgs[s][0]].stddev; // Takes the first RG of the sample as reference.
                const auto currentClusterIt = clusterMaps[s].clusterMap.find(*chromIt);
                if (currentClusterIt != clusterMaps[s].clusterMap.end())
                {
                    // TODO: Improve by considering the individual positional clusters, not only the max cluster.
                    // Strings for holding the pos of supporting reads by orientation and for both loca/remote chrom.
                    String<unsigned> localFRpositions;
                    String<unsigned> localRFpositions;
                    String<unsigned> localFFpositions;
                    String<unsigned> localRRpositions;
                    String<unsigned> remoteFRpositions;
                    String<unsigned> remoteRFpositions;
                    String<unsigned> remoteFFpositions;
                    String<unsigned> remoteRRpositions;

                    call.perSampleTranslocSupport[s].i[0] += getMaxPosCluster(localFRpositions,
                                                                              remoteFRpositions,
                                                                              currentClusterIt,
                                                                              stddev,
                                                                              Orientation::FR,
                                                                              3);
                    call.perSampleTranslocSupport[s].i[1] += getMaxPosCluster(localRFpositions,
                                                                              remoteRFpositions,
                                                                              currentClusterIt,
                                                                              stddev,
                                                                              Orientation::RF,
                                                                              3);
                    call.perSampleTranslocSupport[s].i[2] += getMaxPosCluster(localFFpositions,
                                                                              remoteFFpositions,
                                                                              currentClusterIt,
                                                                              stddev,
                                                                              Orientation::FF,
                                                                              3);
                    call.perSampleTranslocSupport[s].i[3] += getMaxPosCluster(localRRpositions,
                                                                              remoteRRpositions,
                                                                              currentClusterIt,
                                                                              stddev,
                                                                              Orientation::RR,
                                                                              3);
                    call.perSampleTranslocSupport[s].i[4] += sum(sampleWiseOrientations[s]);

                    double p = 0.8;
                    if (!empty(localFRpositions))
                    {
                        call.perSampleFRPositions[s].i1 = getPercentile(localFRpositions, p, Sorted());
                        ++orientationCounter.i[0];
                        ++totalFrRfEvidence;
                    }
                    if (!empty(remoteFRpositions))
                    {
                        call.perSampleFRPositions[s].i2 = getPercentile(remoteFRpositions, 1 - p, Sorted());
                        ++totalFrRfEvidence;
                    }
                    if (!empty(localRFpositions))
                    {
                        ++orientationCounter.i[1];
                        call.perSampleRFPositions[s].i1 = getPercentile(localRFpositions, 1 - p, Sorted());
                        ++totalFrRfEvidence;
                    }
                    if (!empty(remoteRFpositions))
                    {
                        call.perSampleRFPositions[s].i2 = getPercentile(remoteRFpositions, p, Sorted());
                        ++totalFrRfEvidence;
                    }
                    if (!empty(localFFpositions))
                    {
                        ++orientationCounter.i[2];
                        call.perSampleFFPositions[s].i1 = getPercentile(localFFpositions, p, Sorted());
                        ++totalFfRrEvidence;
                    }
                    if (!empty(remoteFFpositions))
                    {
                        call.perSampleFFPositions[s].i2 = getPercentile(remoteFFpositions, p, Sorted());
                        ++totalFfRrEvidence;
                    }
                    if (!empty(localRRpositions))
                    {
                        ++orientationCounter.i[3];
                        call.perSampleRRPositions[s].i1 = getPercentile(localRRpositions, 1 - p, Sorted());
                        ++totalFfRrEvidence;
                    }
                    if (!empty(remoteRRpositions))
                    {
                        call.perSampleRRPositions[s].i2 = getPercentile(remoteRRpositions, 1 - p, Sorted());
                        ++totalFfRrEvidence;
                    }
                    // std::cout << "Adding " << call.perSampleTranslocSupport[s].i[0]
                    //           << " reads to per sample FR support for sample " << s << "." << std::endl;
                    // std::cout << "Adding " << call.perSampleTranslocSupport[s].i[1]
                    //           << " reads to per sample RF support for sample " << s << "." << std::endl;
                    // std::cout << "Adding " << call.perSampleTranslocSupport[s].i[2]
                    //           << " reads to per sample FF support for sample " << s << "." << std::endl;
                    // std::cout << "Adding " << call.perSampleTranslocSupport[s].i[3]
                    //           << " reads to per sample RR support for sample " << s << "." << std::endl;
                    // std::cout << "Adding " << call.perSampleTranslocSupport[s].i[4]
                    //           << " reads to per sample non support for sample " << s << "." << std::endl;
                    // std::cout << "FR(local, remote):" << call.perSampleFRPositions[s] << ", RF(local, remote):"
                    //           << call.perSampleRFPositions[s] << "FF(local, remote):" << call.perSampleFFPositions[s]
                    //           << ", RR(local, remote):" << call.perSampleRRPositions[s] << std::endl;
                }
            }
        }
        assignTranslocOrientation(call, orientationCounter);
        Pair<unsigned> localRemotePos;
        if (totalFrRfEvidence >= totalFfRrEvidence)
            localRemotePos = getLocalAndRemotePos(call.perSampleFRPositions, call.perSampleRFPositions);
        else
            localRemotePos = getLocalAndRemotePos(call.perSampleFFPositions, call.perSampleRRPositions);

        call.position = localRemotePos.i1;
        call.matePosition = localRemotePos.i2;
        call.totalSupport = 1;      // Placeholder value s.t. the call is not marked as invalid.
        //std::cout << "Pos(locals/remote):\t" << call.position << " + " <<  call.matePosition << std::endl;
        appendValue(calls, call);
    }
    return true;
}
// -----------------------------------------------------------------------------
// Function compute_sample_translocation_likelihoods()
// -----------------------------------------------------------------------------
// Calculate the likelihoods for the three genotypes fot the given breakpoint
// TODO: Simple prototype. Refine further
inline void compute_sample_translocation_likelihoods(Triple<double> & sampleGtLogs,
                                                     const unsigned translocationReadCount,
                                                     const unsigned nonTranslocationReadCount,
                                                     //const unsigned breakpoint,
                                                     const double pseudo = 0.001)
{
    sampleGtLogs.i1 = log10(1.0 - pseudo) * nonTranslocationReadCount + log10(pseudo) * translocationReadCount;
    sampleGtLogs.i2 = -log10(2) * (nonTranslocationReadCount + translocationReadCount);
    sampleGtLogs.i3 = log10(pseudo) * nonTranslocationReadCount + log10(1.0 - pseudo) * translocationReadCount;
    //std::cout << sampleGtLogs.i1 << ", " << sampleGtLogs.i2 << ", " << sampleGtLogs.i3 << std::endl;
}
// -----------------------------------------------------------------------------
// Function translocation_likelihood_ratio()
// -----------------------------------------------------------------------------
// Calculate and return the translocation likekelihood ratio from the genotype LOG likelihoods of a single sample.
inline double translocation_likelihood_ratio(const Triple<double> & sampleGtLogs)
{
    double trans = log10(pow(10, sampleGtLogs.i1) + pow(10, sampleGtLogs.i2) + pow(10, sampleGtLogs.i3));
    double no_trans = sampleGtLogs.i1;
    return trans - no_trans;
}
// -----------------------------------------------------------------------------
// Function checkTranslocationClusterThresholds()
// -----------------------------------------------------------------------------
// Check all clusters of the given ClusterMap against the thresholds and set their filter value accordingly.
// Return true if any of the cluster passes the filters, false otherwise.
inline bool checkTranslocationClusterThresholds(TranslocationClusterMap & clusteredRecords,
                                                const unsigned normalCount,
                                                const unsigned absThreshold,
                                                const unsigned relativeThreshold)
{
    bool anyPass = false;
    for (auto it = clusteredRecords.begin(); it != clusteredRecords.end(); ++it)
    {
        unsigned maxCount = std::max(std::max(getClusterSize(it, Orientation::FR),
                                              getClusterSize(it, Orientation::RF)),
                                     std::max(getClusterSize(it, Orientation::FF),
                                              getClusterSize(it, Orientation::RR)));
        SEQAN_ASSERT_GT(maxCount, 0u);
        if (maxCount < absThreshold)
            continue;
        else if (static_cast<double>(maxCount) / (normalCount + getClusterSize(it)) < relativeThreshold)
            continue;
        else
        {
            setClusterFilter(it, true);
            anyPass = true;
        }
    }
    return anyPass;
}
// -----------------------------------------------------------------------------
// Function checkTranslocationRecords()
// -----------------------------------------------------------------------------
// Return true if the reads in the curent window of THIS ONE SAMPLE pass the conditions for translocation counts.
// Return false otherwise.
// Also populate the translocation cluster map for the sample.
inline bool checkTranslocationRecords(const TranslocationProfile & translocProfiles,
                                      TranslocationClusterMap & clusteredRecords,
                                      const unsigned nonTranslocationReadCount,
                                      const TRGs & rgs,
                                      const unsigned sampleIndex,
                                      const unsigned absThreshold = 2,
                                      const double relativeThreshold = 0.1)
{
    unsigned totalTranslocationCount = translocProfiles.getTotalTranslocationReadCounts(rgs[sampleIndex]);
    unsigned totalCount = totalTranslocationCount + nonTranslocationReadCount;
    if ((totalTranslocationCount < absThreshold) ||
        (static_cast<double>(totalTranslocationCount) / totalCount < relativeThreshold))
    {
        return false;
    }
    else    // Cluster by chromosome check counts per orientation.
    {
        clusterTranslocationRecords(clusteredRecords, translocProfiles, rgs[sampleIndex]);
        return checkTranslocationClusterThresholds(clusteredRecords,
                                                   nonTranslocationReadCount,
                                                   absThreshold,
                                                   relativeThreshold);
    }
}
// -----------------------------------------------------------------------------
// Function getAnyPasses()
// -----------------------------------------------------------------------------
// Return a string of bools, each indicating if the respective counts of RF/FF/RR read pairs of ANY SAMPLE
// passed the thresholds.
inline String<bool> getAnyPasses(const String<String<bool> > & sampleWisePasses)
{
    unsigned n = length(sampleWisePasses);
    String<bool> passes;
    resize(passes, 3, false, Exact());
    for (unsigned i = 0; i < 3; ++i)
    {
        for (unsigned s = 0; s < n; ++s)
        {
            if (sampleWisePasses[s][i])
            {
                passes[i] = true;
                break;
            }
        }
    }
    return passes;
}
// -----------------------------------------------------------------------------
// Function getAnyTranslocationPass()
// -----------------------------------------------------------------------------
// Return a bool indication if the any sample in the string is indicated to have a pass.
inline bool getAnyTranslocationPass(const String<bool> & sampleWiseTranslocationPasses)
{
    for (auto it = begin(sampleWiseTranslocationPasses); it != end(sampleWiseTranslocationPasses); ++it)
        if (*it)
            return true;

    return false;
}
// -----------------------------------------------------------------------------
// Function processRFPasses()
// -----------------------------------------------------------------------------
// Return false if the RF oriented read pairs did not pass the filters (indicated by bool frPass).
// Fill all information into call otherwise and return true.
inline bool processRFPasses(JunctionCall & call,
                            const ChromosomeProfile & chromosomeProfiles,
                            const TRGs & rgs,
                            const bool & frPass,
                            const String<String<bool> > & sampleWisePasses,
                            const unsigned & rfCount,
                            const String<OrientationCounts> & sampleWiseOrientations)
{
    if (!frPass)
        return false;

    const unsigned sampleNum = length(rgs);
    call.orientation = Orientation::RF;
    call.svtype=SVType::DUP;
    call.totalSupport = rfCount;
    resize(call.perSampleSupport, sampleNum, 0, Exact());
    Pair<String<unsigned> > perSampleFirstLast;
    for (unsigned s = 0; s < sampleNum; ++s)
    {
        if (sampleWisePasses[s][0])
        {
            call.perSampleSupport[s] = sampleWiseOrientations[s].rfCount;
            Pair<unsigned> fl = chromosomeProfiles.getActiveReadsFirstLast(rgs[s], call.orientation);
            if (fl.i1 != 0)
                appendValue(perSampleFirstLast.i1, fl.i1);
            if (fl.i2 != 0)
                appendValue(perSampleFirstLast.i2, fl.i2);
        }
    }
    double p = 0.2;
    call.position = length(perSampleFirstLast.i1) > 0 ? getPercentile(perSampleFirstLast.i1, p) : 0;
    call.matePosition = length(perSampleFirstLast.i2) > 0 ? getPercentile(perSampleFirstLast.i2, 1.0 - p ) : 0;
    return true;
}
// -----------------------------------------------------------------------------
// Function processFFPasses()
// -----------------------------------------------------------------------------
// Return false if the FF oriented read pairs did not pass the filters (indicated by bool ffPass).
// Fill all information into call otherwise and return true.
inline bool processFFPasses(JunctionCall & call,
                            const ChromosomeProfile & chromosomeProfiles,
                            const TRGs & rgs,
                            const bool & ffPass,
                            const String<String<bool> > & sampleWisePasses,
                            const unsigned & ffCount,
                            const String<OrientationCounts> & sampleWiseOrientations)
{
    if (!ffPass)
        return false;

    const unsigned sampleNum = length(rgs);
    call.orientation = Orientation::FF;
    call.svtype=SVType::INV;
    call.totalSupport = ffCount;
    resize(call.perSampleSupport, sampleNum, 0, Exact());
    Pair<String<unsigned> > perSampleFirstLast;
    for (unsigned s = 0; s < sampleNum; ++s)
    {
        if (sampleWisePasses[s][1])
        {
            call.perSampleSupport[s] = sampleWiseOrientations[s].ffCount;
            Pair<unsigned> fl = chromosomeProfiles.getActiveReadsMaxFirstLast(rgs[s], call.orientation);
            if (fl.i1 != 0)
                appendValue(perSampleFirstLast.i1, fl.i1);
            if (fl.i2 != 0)
                appendValue(perSampleFirstLast.i2, fl.i2);
        }
    }
    double p = 0.8;
    call.position = length(perSampleFirstLast.i1) > 0 ? getPercentile(perSampleFirstLast.i1, p) : 0;
    call.matePosition = length(perSampleFirstLast.i2) > 0 ? getPercentile(perSampleFirstLast.i2, p ) : 0;
    return true;
}
// -----------------------------------------------------------------------------
// Function processRRPasses()
// -----------------------------------------------------------------------------
// Return false if the RR oriented read pairs did not pass the filters (indicated by bool rrPass).
// Fill all information into call otherwise and return true.
inline bool processRRPasses(JunctionCall & call,
                            const ChromosomeProfile & chromosomeProfiles,
                            const TRGs & rgs,
                            const bool & rrPass,
                            const String<String<bool> > & sampleWisePasses,
                            const unsigned & rrCount,
                            const String<OrientationCounts> & sampleWiseOrientations)
{
    if (!rrPass)
        return false;

    const unsigned sampleNum = length(rgs);
    call.orientation = Orientation::RR;
    call.svtype=SVType::INV;
    call.totalSupport = rrCount;
    resize(call.perSampleSupport, sampleNum, 0, Exact());
    Pair<String<unsigned> > perSampleFirstLast;
    for (unsigned s = 0; s < sampleNum; ++s)
    {
        if (sampleWisePasses[s][2])
        {
            call.perSampleSupport[s] = sampleWiseOrientations[s].rrCount;
            Pair<unsigned> fl = chromosomeProfiles.getActiveReadsMinFirstLast(rgs[s], call.orientation);
            if (fl.i1 != 0)
                appendValue(perSampleFirstLast.i1, fl.i1);
            if (fl.i2 != 0)
                appendValue(perSampleFirstLast.i2, fl.i2);
        }
    }
    double p = 0.2;
    call.position = length(perSampleFirstLast.i1) > 0 ? getPercentile(perSampleFirstLast.i1, p) : 0;
    call.matePosition = length(perSampleFirstLast.i2) > 0 ? getPercentile(perSampleFirstLast.i2, p ) : 0;
    return true;
}
// -----------------------------------------------------------------------------
// Function assignCoverageChanges()
// -----------------------------------------------------------------------------
// Take ChromosomeProfile::previousActiveLoad, calculate the current active coverage and set the
// coverage change values of the call accordingly.
// Update ChromosomeProfile::previousActiveLoad with the current active coverage when done.
inline void assignCoverageChanges(String<JunctionCall> & rfCalls,
                                  String<JunctionCall> & ffCalls,
                                  String<JunctionCall> & rrCalls,
                                  ChromosomeProfile & chromosomeProfiles,
                                  const TRGs & rgs,
                                  const bool & callMade)
{
    if (callMade)
    {
        if (!empty(rfCalls))
            if (back(rfCalls).windowPosition == chromosomeProfiles.currentPos)
                chromosomeProfiles.assignAndupdateActiveLoad(back(rfCalls), rgs);
        if (!empty(ffCalls))
            if (back(ffCalls).windowPosition == chromosomeProfiles.currentPos)
                chromosomeProfiles.assignAndupdateActiveLoad(back(ffCalls), rgs);
        if (!empty(rrCalls))
            if (back(rrCalls).windowPosition == chromosomeProfiles.currentPos)
                chromosomeProfiles.assignAndupdateActiveLoad(back(rrCalls), rgs);
    }
    else
    {
        chromosomeProfiles.updateActiveLoad(rgs);   // Update ActiveLoad, but don't assign changes to a call.
    }
    //std::cout << chromosomeProfiles.currentPos << "\t" << chromosomeProfiles.previousActiveLoad[0] << std::endl;
}
// -----------------------------------------------------------------------------
// Function compute_sample_duplication_likelihoods()
// -----------------------------------------------------------------------------
// Calculate the log likelihood for all three genotypes with the given number of supporting an non support read pairs.
inline void compute_sample_duplication_likelihoods(Triple<double> & sampleGtLogs,
                                                   const ChromosomeProfile & chromosomeProfiles,
                                                   const String<unsigned> & sampleRGs,
                                                   const String<Histogram> & hists,
                                                   const unsigned breakpointOne,
                                                   const unsigned breakpointTwo)
{
    for (const unsigned & rg : sampleRGs)
    {
        const Histogram & hist = hists[rg];
        ChromosomeProfile::TActiveSet::const_iterator it(chromosomeProfiles.activeReads[rg].begin());
        ChromosomeProfile::TActiveSet::const_iterator itEnd(chromosomeProfiles.activeReads[rg].end());
        // std::cout << "====Examining new read group====\nSV: " << breakpointOne  << "-" << breakpointTwo
        //            << "(" << breakpointTwo - breakpointOne << ")" << std::endl;
        // unsigned refSup = 0;
        // unsigned dupSup= 0;
        // unsigned bothSup = 0;

        for (;it != itEnd; ++it)
        {
            // std::cout << "Inspecting record " << chromosomeProfiles.getSingleStartPos(rg, it) << "-"
            //           << chromosomeProfiles.getSingleEndPos(rg, it) << ":"
            //           << chromosomeProfiles.getSingleOrientation(rg, it) << std::endl;
            int currentDeviation = chromosomeProfiles.getSingleDeviation(rg, it);
            // std::cout << "deviation = " << currentDeviation <<  std::endl;
            long double dupLikelihood = I(hist, currentDeviation - (breakpointTwo - breakpointOne));
            long double refLikelihood = I(hist, currentDeviation);
            // if (dupLikelihood > refLikelihood)
            //     ++dupSup ;
            // else if (refLikelihood > dupLikelihood)
            //     ++refSup;
            // else
            //     ++bothSup;
            // std::cout << "Lref: " << refLikelihood << "\tLDup: " << dupLikelihood << std::endl;

            // std::cout << "G0 = " << sampleGtLogs.i1 << " + " << log10(refLikelihood) << " = " ;
            sampleGtLogs.i1 += log10(refLikelihood);
            // std::cout << sampleGtLogs.i1 << std::endl;
            // std::cout << "G1 = " << sampleGtLogs.i2 << " + " << log10(2.0/3.0 * refLikelihood + 1.0/3.0 * dupLikelihood) << " = ";
            sampleGtLogs.i2 += log10(0.5 * refLikelihood + 0.5 * dupLikelihood);
            // std::cout << sampleGtLogs.i2 << std::endl;
            // std::cout << "G2 = " << sampleGtLogs.i3 << " + " << log10(0.5 * refLikelihood + 0.5 * dupLikelihood) << " = ";
            sampleGtLogs.i3 += log10(1.0/3.0 * refLikelihood + 2.0/3.0 * dupLikelihood);
            // std::cout << sampleGtLogs.i3 << std::endl;
        }
        // std::cout << "Final G0;G1,G2 = " << sampleGtLogs.i1 << ";" << sampleGtLogs.i2 << ";" << sampleGtLogs.i3 << std::endl;
        // std::cout << "RefSup = " << refSup << "\tDupSup = " << dupSup << "\tAmbigous = " << bothSup << std::endl;
    }
}
// -----------------------------------------------------------------------------
// Function pairSpansSingleBreakpoint()
// -----------------------------------------------------------------------------
// Return true if a read pair spans exactly one breakpoint, i.e. one end lies within the variant, the other one outside.
inline bool pairSpansSingleBreakpoint(const unsigned pairStart,
                                      const unsigned pairEnd,
                                      const unsigned breakpointOne,
                                      const unsigned breakpointTwo,
                                      const unsigned readLength,
                                      const Orientation orientation)
{
    unsigned left = pairStart;
    unsigned right = pairEnd;
    if (orientation == Orientation::FR)
    {
        left -= readLength / 2;
        right += readLength / 2;
    }
    else if (orientation == Orientation::FF)
    {
        left -= readLength / 2;
        right -= readLength / 2;
    }
    else if (orientation == Orientation::RR)
    {   left += readLength / 2;
        right += readLength / 2;
    }
    if (left <= breakpointOne)
        return right <= breakpointTwo;
    else if (left <= breakpointTwo)
        return right > breakpointTwo;
    else
        return false;
}
// -----------------------------------------------------------------------------
// Function compute_sample_inv_likelihoods()
// -----------------------------------------------------------------------------
// Calculate the log likelihood for all three genotypes with the given number of supporting an non support read pairs.
inline void compute_sample_inv_likelihoods(Triple<double> & sampleGtLogs,
                                           const ChromosomeProfile & chromosomeProfiles,
                                           const String<unsigned> & sampleRGs,
                                           const String<Histogram> & hists,
                                           const unsigned breakpointOne,
                                           const unsigned breakpointTwo,
                                           const double pseudo = 0.00001)
{
    for (const unsigned & rg : sampleRGs)
    {
        const Histogram & hist = hists[rg];
        ChromosomeProfile::TActiveSet::const_iterator it(chromosomeProfiles.activeReads[rg].begin());
        ChromosomeProfile::TActiveSet::const_iterator itEnd(chromosomeProfiles.activeReads[rg].end());
        // std::cout << "====Examining new read group====\nSV: " << breakpointOne  << "-" << breakpointTwo
        //           << "(" << breakpointTwo - breakpointOne << ")" << std::endl;
        for (;it != itEnd; ++it)
        {
            Orientation orientation = chromosomeProfiles.getSingleOrientation(rg, it);
            int left = chromosomeProfiles.getSingleStartPos(rg, it);
            int right = chromosomeProfiles.getSingleEndPos(rg, it);
            if (pairSpansSingleBreakpoint(left,
                                          right,
                                          breakpointOne,
                                          breakpointTwo,
                                          hists[rg].readLength,
                                          orientation))
            {
                int a = std::abs(left - static_cast<int>(breakpointOne));
                int b = std::abs(right - static_cast<int>(breakpointTwo));
                int uncorrectedDeviation = chromosomeProfiles.getSingleDeviation(rg, it);
                int correctedDeviation;
                if (orientation == Orientation::FF || orientation == Orientation::RR)
                    correctedDeviation = a + b - hists[rg].median3PrimeDist - chromosomeProfiles.getSingleClipping(rg, it);
                else
                    correctedDeviation = uncorrectedDeviation;

                // std::cout << "Orientation = " << orientation << std::endl;
                // std::cout << "a, b, clip, median = " << a << ", " << b << ", " << chromosomeProfiles.getSingleClipping(rg, it) << "," <<  hists[rg].median3PrimeDist << std::endl;
                // std::cout << "uncorrectedDeviation = " << uncorrectedDeviation <<  std::endl;
                // std::cout << "correcteddeviation = " << correctedDeviation <<  std::endl;

                // std::cout << "Inspecting record " << chromosomeProfiles.getSingleStartPos(rg, it) << "-"
                //         << chromosomeProfiles.getSingleEndPos(rg, it) << ":"
                //         << orientation << std::endl;

                long double refLikelihood = I(hist, uncorrectedDeviation);
                long double invLikelihood = pseudo;
                 if (orientation == Orientation::FF || orientation == Orientation::RR)
                    invLikelihood = I(hist, correctedDeviation);

                sampleGtLogs.i1 += log10(refLikelihood);
                sampleGtLogs.i2 += log10(refLikelihood + invLikelihood) - log10(2.0);
                sampleGtLogs.i3 += log10(invLikelihood);
                // std::cout << refLikelihood << "\t" << invLikelihood << std::endl;
                // std::cout << "Local G0;G1,G2 = " << log10(refLikelihood) << ";" << log10(refLikelihood + invLikelihood) - log10(2.0) << ";" << log10(invLikelihood) << std::endl;
                // std::cout << "Updt. G0;G1,G2 = " << sampleGtLogs.i1 << ";" << sampleGtLogs.i2 << ";" << sampleGtLogs.i3 << std::endl;

                // if (orientation == Orientation::FR)
                // {
                //     long double refLikelihood = I(hist, currentDeviation);
                //     sampleGtLogs.i1 += log10(refLikelihood) + log10(2.0/3.0);
                //     sampleGtLogs.i2 += log10(refLikelihood) + log10(1.0/3.0);
                //     sampleGtLogs.i3 += log10(pseudo);
                // }
                // else if (orientation == Orientation::FF || orientation == Orientation::RR)
                // {
                //     long double invLikelihood = I(hist, currentDeviation - (breakpointTwo - breakpointOne));
                //     sampleGtLogs.i1 += log10(pseudo);
                //     sampleGtLogs.i2 += log10(invLikelihood) + log10(1.0/3.0);
                //     sampleGtLogs.i3 += log10(invLikelihood) + log10(2.0/3.0);
                // }
                // else
                // {
                //     long double refLikelihood = I(hist, currentDeviation);
                //     sampleGtLogs.i1 += log10(refLikelihood) + log10(pseudo);
                //     sampleGtLogs.i2 += log10(refLikelihood) + log10(pseudo);
                //     sampleGtLogs.i3 += log10(refLikelihood) + log10(pseudo);
                // }
                // if (orientation == Orientation::FR)
                // {
                //     sampleGtLogs.i1 += log10(2.0/3.0);
                //     sampleGtLogs.i2 += log10(1.0/3.0);
                //     sampleGtLogs.i3 += log10(pseudo);
                // }
                // else if (orientation == Orientation::FF || orientation == Orientation::RR)
                // {
                //     sampleGtLogs.i1 += log10(pseudo);
                //     sampleGtLogs.i2 += log10(1.0/3.0);
                //     sampleGtLogs.i3 += log10(2.0/3.0);
                // }
                // else
                // {
                //     sampleGtLogs.i1 += log10(pseudo);
                //     sampleGtLogs.i2 += log10(pseudo);
                //     sampleGtLogs.i3 += log10(pseudo);
                // }
            }
            // else
            // {
            //             std::cout << "Record " << chromosomeProfiles.getSingleStartPos(rg, it) << "-"
            //             << chromosomeProfiles.getSingleEndPos(rg, it) << ":"
            //             << chromosomeProfiles.getSingleOrientation(rg, it)
            //             << " is not overlapping the breakpoints and will be ignored." << std::endl;
            // }
        }
        // std::cout << "Final G0;G1,G2 = " << sampleGtLogs.i1 << ";" << sampleGtLogs.i2 << ";" << sampleGtLogs.i3 << std::endl;
    }
}
// -----------------------------------------------------------------------------
// Function phredTransform()
// -----------------------------------------------------------------------------
// Perform PHRED scaling on a triple of log-likelihoods.
inline void phredTransform(Triple<double> & gtLogs)
{
    double g0 = gtLogs.i1;
    double g1 = gtLogs.i2;
    double g2 = gtLogs.i3;
    const double gTot = log10(exp(g0) + exp(g1) + exp(g2)); // g0-g2 where calculated using ln.
    g0 = -10 * (g0 - gTot);
    g1 = -10 * (g1 - gTot);
    g2 = -10 * (g2 - gTot);
    const double minPhred = std::min(std::min(g0, g1), g2);
    gtLogs.i1 = round(g0 - minPhred);
    gtLogs.i2 = round(g1 - minPhred);
    gtLogs.i3 = round(g2 - minPhred);
}
// -----------------------------------------------------------------------------
// Function checkPerSampleTranslocationFilters()
// -----------------------------------------------------------------------------
// Check if all additional filters for translocations are beeing passed for a single sample.
// Return false otherwise.
inline bool checkPerSampleTranslocationFilters(const TranslocationClusterMap & translocClusterMap,
                                               const unsigned posClustersFR,
                                               const unsigned posClustersRF,
                                               const unsigned posClustersFF,
                                               const unsigned posClustersRR,
                                               const unsigned maxPasses = 1u,
                                               const unsigned maxClusters = 2u,
                                               const unsigned maxPosClusters = 2u)
{
    if (translocClusterMap.getClusterPasses() > maxPasses)
    {
        std::cout << "Ignoring window with too many (" << translocClusterMap.getClusterPasses()
                  << ") passing translocation clusters" << std::endl;
        return false;
    }
    else if (translocClusterMap.size() > maxClusters)
    {
        std::cout << "Ignoring window with too many ("<< translocClusterMap.size()
                  << ") translocation clusters" << std::endl;
        return false;
    }
    else if (posClustersFR > maxPosClusters)
    {
        std::cout << "Ignoring window with too many (" << posClustersFR << ") FR positional clusters." << std::endl;
        return false;
    }
    else if (posClustersRF > maxPosClusters)
    {
        std::cout << "Ignoring window with too many (" << posClustersRF << ") RF positional clusters." << std::endl;
        return false;
    }
    else if (posClustersFF > maxPosClusters)
    {
        std::cout << "Ignoring window with too many (" << posClustersFF << ") FF positional clusters." << std::endl;
        return false;
    }
    else if (posClustersRR > maxPosClusters)
    {
        std::cout << "Ignoring window with too many (" << posClustersRR << ") RR positional clusters." << std::endl;
        return false;
    }
    else if ((posClustersRF >= 1 || posClustersFR >= 1) && (posClustersFF >= 1 || posClustersRR >= 1))
    {
        std::cout << "Ignoring window with both (FR or RF) and (FF or RR) clusters." << std::endl;
        return false;
    }
    else
    {
        return true;
    }
}
// -----------------------------------------------------------------------------
// Function checkAndMarkBadTranslocationSample()
// -----------------------------------------------------------------------------
// Check the per sample translocation filters and mark failing samples as bad. 
inline void checkAndMarkBadTranslocationSample(String<TranslocationClusterMap> & translocClusterMaps,
                                               const String<bool> & sampleWiseTranslocationPasses,
                                               const String<Histogram> & histograms,
                                               const TRGs & rgs,
                                               const unsigned s)         // Sample index/number
{
    if (sampleWiseTranslocationPasses[s])
    {
        double stddev = histograms[rgs[s][0]].stddev; // Takes the first RG of the sample as reference
        unsigned posClustersFR = 0;
        unsigned posClustersRF = 0;
        unsigned posClustersFF = 0;
        unsigned posClustersRR = 0;
        for (auto it = translocClusterMaps[s].begin(); it != translocClusterMaps[s].end(); ++it)
        {
            posClustersFR += getNumPositionClusters(it, stddev, Orientation::FR, 3);
            posClustersRF += getNumPositionClusters(it, stddev, Orientation::RF ,3);
            posClustersFR += getNumPositionClusters(it, stddev, Orientation::FF, 3);
            posClustersRF += getNumPositionClusters(it, stddev, Orientation::RR ,3);
        }
        if (checkPerSampleTranslocationFilters(translocClusterMaps[s],
                                               posClustersFR,
                                               posClustersRF,
                                               posClustersFF,
                                               posClustersRR))
        {
            //For testing: Output of passing clusteers with preliminary likelihoods
            // for (auto it = translocClusterMaps[s].begin(); it != translocClusterMaps[s].end(); ++it)
            // {
            //     Triple<double> translocationGenotypeLikelihoods; // Just for testing
            //     unsigned translocationReadCount = getClusterSize(it, Orientation::FR) +
            //                                       getClusterSize(it, Orientation::RF);
            //     compute_sample_translocation_likelihoods(translocationGenotypeLikelihoods,
            //                                             translocationReadCount,
            //                                             sum(sampleWiseOrientations[s]));
            //     if (translocationGenotypeLikelihoods.i2 > translocationGenotypeLikelihoods.i1 ||
            //         translocationGenotypeLikelihoods.i3 > translocationGenotypeLikelihoods.i1)
            //     {
            //         std::cout << "Potential translocation junction in sample " << s << " in window "
            //                 << chromosomeProfiles.chrom << "\t" << chromosomeProfiles.currentPos << "\t"
            //                 << chromosomeProfiles.currentPos + 30 << " Sup = "
            //                 << translocClusterMaps[s].begin()->second.i1.size() << "/"
            //                 << sum(sampleWiseOrientations[s])
            //                 << " with " << translocClusterMaps[s].getClusterPasses()  << " passing cluster(s) ("
            //                 <<  translocClusterMaps[s].size()  << " in total)."
            //                 << std::endl;
            //         std::cout << translocationGenotypeLikelihoods << std::endl;
            //         std::cout << translocation_likelihood_ratio(translocationGenotypeLikelihoods) << std::endl;
            //     }
            //     else
            //     {
            //         //std::cout << "Ignoring window with L(G0) > L(G[1/2])" << std::endl;
            //     }
            // }
        }
        else    // Mark the sample as bad by setting all of it's clusters to non-passing
        {
            translocClusterMaps[s].markAsBadSample();
        }
    }
}
// -----------------------------------------------------------------------------
// Function detect_junction_window()
// -----------------------------------------------------------------------------
// Inspect the read pairs of the current window for potential novel junctions.
// Return true if at least one has been found, false otherwise.
inline bool detect_junction_window(String<JunctionCall> & rfCalls,
                                   String<JunctionCall> & ffCalls,
                                   String<JunctionCall> & rrCalls,
                                   String<JunctionCall> & translocCalls,
                                   const ChromosomeProfile & chromosomeProfiles,
                                   const TranslocationProfile & translocProfiles,
                                   const TRGs & rgs,
                                   PopDelCallParameters & params)
{
    unsigned sampleCount = length(rgs);
    String<OrientationCounts> sampleWiseOrientations;
    resize(sampleWiseOrientations, sampleCount, OrientationCounts(0, 0, 0, 0), Exact());
    String<String<bool> > sampleWisePasses;
    resize(sampleWisePasses, sampleCount, Exact());
    // String<bool> sampleWiseTranslocationPasses;
    // resize(sampleWiseTranslocationPasses, sampleCount, Exact());
    OrientationCounts totalOrientationCounts(0, 0, 0, 0);
    // unsigned totalTranslocationCount = 0;
    // String<TranslocationClusterMap> translocClusterMaps;        // For each sample a map<chrom, records>
    // resize(translocClusterMaps, sampleCount, Exact());
    for (unsigned s = 0; s < sampleCount; ++s)          // For each sample
    {
        sampleWiseOrientations[s] += chromosomeProfiles.getActiveReadsOrientations(rgs[s]);
        totalOrientationCounts += sampleWiseOrientations[s];
        // totalTranslocationCount += translocClusterMaps[s].size();
        sampleWisePasses[s] = checkOrientationCounts(sampleWiseOrientations[s],
                                                     params.orientationThreshold,
                                                     params.relOrientationThreshold);
        // sampleWiseTranslocationPasses[s] = checkTranslocationRecords(translocProfiles,
        //                                                              translocClusterMaps[s],
        //                                                              sum(sampleWiseOrientations[s]),
        //                                                              rgs,
        //                                                              s,
        //                                                              params.translocationThreshold,
        //                                                              params.relTranslocationThreshold);
        // checkAndMarkBadTranslocationSample(translocClusterMaps,
        //                                    sampleWiseTranslocationPasses,
        //                                    params.histograms,
        //                                    rgs,
        //                                    s);
    }
    String<bool> totalPasses = getAnyPasses(sampleWisePasses);   // TODO: Make this a tuple instead of a string.
    // bool anyTranslocationPass = getAnyTranslocationPass(sampleWiseTranslocationPasses);
    // bool translocProcess = processTranslocationClusters(translocCalls,
    //                                                     translocClusterMaps,
    //                                                     sampleWiseOrientations,
    //                                                     rgs,
    //                                                     sampleWiseTranslocationPasses,
    //                                                     anyTranslocationPass,
    //                                                     chromosomeProfiles.currentPos,
    //                                                     params);

    JunctionCall rfJunction(chromosomeProfiles.currentPos);                                   // TODO: Multiple inits for overlapping SVs of one type.
    bool rfProcess = processRFPasses(rfJunction,
                                     chromosomeProfiles,
                                     rgs,
                                     totalPasses[0],
                                     sampleWisePasses,
                                     totalOrientationCounts.rfCount,
                                     sampleWiseOrientations);
    if (rfProcess)
    {
        bool call = false;
        resize(rfJunction.gtLikelihoods, sampleCount, Exact());
        for (unsigned s = 0; s < sampleCount; ++s)
        {
            compute_sample_duplication_likelihoods(rfJunction.gtLikelihoods[s],
                                                   chromosomeProfiles,
                                                   rgs[s],
                                                   params.histograms,
                                                   rfJunction.position,
                                                   rfJunction.matePosition);
            phredTransform(rfJunction.gtLikelihoods[s]);
            call |= (rfJunction.gtLikelihoods[s].i2 < rfJunction.gtLikelihoods[s].i1) ||
                    (rfJunction.gtLikelihoods[s].i3 < rfJunction.gtLikelihoods[s].i1);
        }
        if (call)
            appendValue(rfCalls, rfJunction);
    }

    JunctionCall ffJunction(chromosomeProfiles.currentPos);
    bool ffProcess = processFFPasses(ffJunction,
                                     chromosomeProfiles,
                                     rgs,
                                     totalPasses[1],
                                     sampleWisePasses,
                                     totalOrientationCounts.ffCount,
                                     sampleWiseOrientations);
    if (ffProcess)
    {
        bool call = false;
        resize(ffJunction.gtLikelihoods, sampleCount, Exact());
        for (unsigned s = 0; s < sampleCount; ++s)
        {
            compute_sample_inv_likelihoods(ffJunction.gtLikelihoods[s],
                                           chromosomeProfiles,
                                           rgs[s],
                                           params.histograms,
                                           ffJunction.position,
                                           ffJunction.matePosition);
            phredTransform(ffJunction.gtLikelihoods[s]);
            call |= (ffJunction.gtLikelihoods[s].i2 < ffJunction.gtLikelihoods[s].i1) ||
                    (ffJunction.gtLikelihoods[s].i3 < ffJunction.gtLikelihoods[s].i1);
        }
        if (call)
            appendValue(ffCalls, ffJunction);
    }
    JunctionCall rrJunction(chromosomeProfiles.currentPos);
    bool rrProcess = processRRPasses(rrJunction,
                                     chromosomeProfiles,
                                     rgs,
                                     totalPasses[2],
                                     sampleWisePasses,
                                     totalOrientationCounts.rrCount,
                                     sampleWiseOrientations);
    if (rrProcess)
    {
        bool call = false;
        resize(rrJunction.gtLikelihoods, sampleCount, Exact());
        for (unsigned s = 0; s < sampleCount; ++s)
        {
            compute_sample_inv_likelihoods(rrJunction.gtLikelihoods[s], // Should yield the same results as for FF
                                           chromosomeProfiles,
                                           rgs[s],
                                           params.histograms,
                                           rrJunction.position,
                                           rrJunction.matePosition);
            phredTransform(rrJunction.gtLikelihoods[s]);
            call |= (rrJunction.gtLikelihoods[s].i2 < rrJunction.gtLikelihoods[s].i1) ||
                    (rrJunction.gtLikelihoods[s].i3 < rrJunction.gtLikelihoods[s].i1);
        }
        if (call)
            appendValue(rrCalls, rrJunction);
    }
    return rfProcess || ffProcess || rrProcess;// || translocProcess;
}
// =======================================================================================
// Function setDupGenotypes()
// =======================================================================================
// Looks at the windows that are within the duplication and estimates the genotype for each sample.
// Resets "genotypes"
// Return true on success, false if no genotype could be set (due to lack of windows)
inline bool setDupGenotypes(const Iterator<String<JunctionCall>, Standard >::Type & start,
                            const Iterator<const String<JunctionCall>, Standard >::Type end,
                            String<Triple<double> > & genotypes)
{
    unsigned genotypeWinCount = 0;
    for (Iterator<String<JunctionCall>, Standard >::Type gtIt = start; gtIt < end; ++gtIt)
    {
        if (gtIt->windowPosition > start->position &&
            gtIt->windowPosition - 30 < start->matePosition)
        {   // Only consider genotype likelihoods of windows within the range of the variant. TODO: Exclude breakpoint windows.
            for (unsigned s = 0; s < length(genotypes); ++s)
            {
                genotypes[s].i1 += gtIt->gtLikelihoods[s].i1;
                genotypes[s].i2 += gtIt->gtLikelihoods[s].i2;
                genotypes[s].i3 += gtIt->gtLikelihoods[s].i3;
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
        double ref = (genotypes[s].i1 - minGt) / genotypeWinCount;
        double het = (genotypes[s].i2 - minGt) / genotypeWinCount;
        double hom = (genotypes[s].i3 - minGt) / genotypeWinCount;
        genotypes[s] = Triple<double> (0.0, 0.0, 0.0);
        start->gtLikelihoods[s].i1 = std::round(ref);
        start->gtLikelihoods[s].i2 = std::round(het);
        start->gtLikelihoods[s].i3 = std::round(hom);
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
// Function setInvGenotypes()
// =======================================================================================
// Looks at the windows that are within the inversion and estimates the genotype for each sample.
// Resets "genotypes"
// Return true on success, false if no genotype could be set (due to lack of windows)
inline bool setInvGenotypes(const Iterator<String<JunctionCall>, Standard >::Type & start,
                            const Iterator<const String<JunctionCall>, Standard >::Type end,
                            String<Triple<double> > & genotypes)
{
    unsigned genotypeWinCount = 0;
    for (Iterator<String<JunctionCall>, Standard >::Type gtIt = start; gtIt < end; ++gtIt)
    {
        if ((gtIt->windowPosition <= start->position &&
             gtIt->windowPosition + 30 > start->position) ||
            (gtIt->windowPosition <= start->matePosition &&
             gtIt->windowPosition + 30 > start->matePosition))
        {   // Only consider breakpoint windows
            for (unsigned s = 0; s < length(genotypes); ++s)
            {
                genotypes[s].i1 += gtIt->gtLikelihoods[s].i1;
                genotypes[s].i2 += gtIt->gtLikelihoods[s].i2;
                genotypes[s].i3 += gtIt->gtLikelihoods[s].i3;
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
        double ref = (genotypes[s].i1 - minGt) / genotypeWinCount;
        double het = (genotypes[s].i2 - minGt) / genotypeWinCount;
        double hom = (genotypes[s].i3 - minGt) / genotypeWinCount;
        genotypes[s] = Triple<double> (0.0, 0.0, 0.0);
        start->gtLikelihoods[s].i1 = std::round(ref);
        start->gtLikelihoods[s].i2 = std::round(het);
        start->gtLikelihoods[s].i3 = std::round(hom);
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
// -----------------------------------------------------------------------------
// Function mergeJunctionCallsRange()
// -----------------------------------------------------------------------------
// Dertermine the position of a junction call from the given strings and assign them the call.
// Reset the containers and counters when done.
inline void mergeJunctionCallsRange(const Iterator<String<JunctionCall> >::Type & anchor,
                                    const Iterator<String<JunctionCall> >::Type & last,
                                    String<unsigned> & positions,
                                    String<unsigned> & matePositions,
                                    String<double> & avgPerSampleSupport,
                                    double & avgTotalSupport,
                                    String<Triple<double> > & genotypes,
                                    unsigned & callCount,
                                    unsigned & winCount)
{
    double p;
    if (anchor->orientation == Orientation::RF)         // Duplication
    {
        p = 0.2;
        anchor->position = getPercentile(positions, p);
        anchor->matePosition = getPercentile(matePositions, 1.0 - p);
        if (setDupGenotypes(anchor, last, genotypes))
        {
            setFreqFromGTs(*anchor);
        }
        else
        {
            invalidateJunctionCall(*anchor);
        }
    }
    else if (anchor->orientation == Orientation::FF)    // Inversion BP1
    {
        p = 0.8;
        anchor->position = getPercentile(positions, p);
        anchor->matePosition = getPercentile(matePositions, p);
        if (setInvGenotypes(anchor, last, genotypes))
        {
            setFreqFromGTs(*anchor);
        }
        else
        {
            invalidateJunctionCall(*anchor);
        }
    }
    else // if (start->orientation == Orientation::RR)  // Inversion BP2
    {
        p = 0.2;
        anchor->position = getPercentile(positions, p);
        anchor->matePosition = getPercentile(matePositions, p);
        if (setInvGenotypes(anchor, last, genotypes))
        {
            setFreqFromGTs(*anchor);
        }
        else
        {
            invalidateJunctionCall(*anchor);
        }
    }
    for(unsigned s = 0; s < length(avgPerSampleSupport); ++s)
        avgPerSampleSupport[s] /= winCount;

    if (isValid(*anchor))
    {
        anchor->perSampleSupport = avgPerSampleSupport;
        avgTotalSupport /= winCount;
        anchor->totalSupport = avgTotalSupport;
        ++callCount;
    }
    clear(positions);
    clear(matePositions);
    for (auto it = begin(avgPerSampleSupport); it != end(avgPerSampleSupport); ++it)
        *it = 0.0;

    winCount = 1;
    avgTotalSupport = 0;
}
// -----------------------------------------------------------------------------
// Operator+= overload
// -----------------------------------------------------------------------------
// Add the values of the second Tuple to the first one
inline Tuple<double, 5> & operator+=(Tuple<double, 5> & l, const Tuple<double, 5> & r)
{
    for (unsigned j = 0; j < 5u; ++j)
        l.i[j] += r.i[j];

    return l;
}

// -----------------------------------------------------------------------------
// Function pairwiseTranslocMerge()
// -----------------------------------------------------------------------------
// Merge the JunctionCall pointed at by source into the one pointed at by target.
// Invalidate the call pointed at by targettarget.
inline void pairwiseTranslocMerge(const Iterator<String<JunctionCall> >::Type & target,
                                  const Iterator<String<JunctionCall> >::Type & source)
{
        std::cout << "Merging: ";
        source->print();
        std::cout << " with ";
        target->print();
        std::cout << std::endl;
        append(target->perSampleFRPositions, source->perSampleFRPositions);
        append(target->perSampleRFPositions, source->perSampleRFPositions);
        append(target->perSampleFFPositions, source->perSampleFFPositions);
        append(target->perSampleRRPositions, source->perSampleRRPositions);
        for (unsigned s = 0; s < length(target->perSampleTranslocSupport); ++s)
        {
            target->perSampleTranslocSupport[s] += source->perSampleTranslocSupport[s];
        }
        invalidateJunctionCall(*source);
}
// -----------------------------------------------------------------------------
// Function assignPerSampleTranslocGTLikelihoods()
// -----------------------------------------------------------------------------
// Assign the GT-likelihoods based on the orientation of the translocation call
inline void assignPerSampleTranslocGTLikelihoods(JunctionCall & c)
{
    resize(c.gtLikelihoods, length(c.perSampleTranslocSupport), Exact());
    if (c.orientation == Orientation::FR || c.orientation == Orientation::RF)
    {
        for (unsigned s = 0; s < length(c.gtLikelihoods); ++s)
        {
            compute_sample_translocation_likelihoods(c.gtLikelihoods[s],
                                                     c.perSampleTranslocSupport[s].i[0] +
                                                        c.perSampleTranslocSupport[s].i[1],
                                                     c.perSampleTranslocSupport[s].i[4]);
            phredTransform(c.gtLikelihoods[s]);
        }
    }
    else
    {
        for (unsigned s = 0; s < length(c.gtLikelihoods); ++s)
        {
            compute_sample_translocation_likelihoods(c.gtLikelihoods[s],
                                                     c.perSampleTranslocSupport[s].i[2] +
                                                        c.perSampleTranslocSupport[s].i[3],
                                                     c.perSampleTranslocSupport[s].i[4]);
            phredTransform(c.gtLikelihoods[s]);
        }
    }
}
// -----------------------------------------------------------------------------
// Function refineMergedTranslocPos()
// -----------------------------------------------------------------------------
inline void refineMergedTranslocPos(JunctionCall & c)
{
    Pair<unsigned> positions;

    std::cout << "Refining: FF:";
    for (auto it = begin(c.perSampleFFPositions); it != end(c.perSampleFFPositions); ++it)
    {
        std::cout << *it << "\t";
    }
    std::cout << "RR: ";
    for (auto it = begin(c.perSampleRRPositions); it != end(c.perSampleRRPositions); ++it)
    {
        std::cout << *it << "\t";
    }
    std::cout << std::endl;

    if (c.orientation == Orientation::FR || c.orientation == Orientation::RF)
    {
        positions = getLocalAndRemotePos(c.perSampleFRPositions, c.perSampleRFPositions);
    }
    else // if(c.orientation == Orientation::FF || c.orientation == Orientation::RR)
    {
        positions = getLocalAndRemotePos(c.perSampleFFPositions, c.perSampleRRPositions);
    }
    c.position = positions.i1;
    c.matePosition = positions.i2;
}
// -----------------------------------------------------------------------------
// Function refineMergedTranslocCall()
// -----------------------------------------------------------------------------
// Rre-calculate the values of a merged translocation call
inline void refineMergedTranslocCall(JunctionCall & c)
{
    assignTranslocOrientation(c);
    // TODO: Filter out calls that have both FR+RF and FF+RR evidence high
    assignPerSampleTranslocGTLikelihoods(c);
    setFreqFromGTs(c);
    refineMergedTranslocPos(c);
}
// -----------------------------------------------------------------------------
// Function refineMergedTranslocCall()
// -----------------------------------------------------------------------------
// Rre-calculate the values of a merged translocation call
// inline void refineMergedTranslocCalls(String<JunctionCall> & junctionCalls & calls)
// {
//     SEQAN_ASSERT(!empty(calls));
//     SEQAN_ASSERT_EQ(calls[0].svtype, SVTYPE::TRL);
//     for (auto it = begin(calls); it != end(calls); ++it)
//         refineMergedTranslocCall(*it);
// }
// -----------------------------------------------------------------------------
// Function mergeJunctionWindows()
// -----------------------------------------------------------------------------
// Merge subsequent junctionCalls into single events.
// Return false if the string of calls is empty, true otherwise.
// NOTE: Should only be applied on string of junctionsCalls of the same orientation!
inline bool mergeJunctionWindows(String<JunctionCall> & junctionCalls, const double & stddev)
{
    if (length(junctionCalls) <= 1u)
        return false;
    const unsigned numRgs = length(junctionCalls[0].perSampleSupport);
    std::sort(begin(junctionCalls), end(junctionCalls));
    String<JunctionCall> tmp;
    Iterator<String<JunctionCall> >::Type currentAnchor = begin(junctionCalls, Standard());
    const Iterator<const String<JunctionCall>, Standard >::Type last = end(junctionCalls) - 1;
    const Iterator<String<JunctionCall> >::Type firstGoodWin = currentAnchor;
    String<unsigned> positions;
    String<unsigned> matePositions;
    String<double> avgPerSampleSupport;
    resize(avgPerSampleSupport, numRgs, 0u, Exact());
    String<Triple<double> > genotypes;
    resize(genotypes, numRgs, Triple<double>(0.0, 0.0, 0.0), Exact());
    append(positions, firstGoodWin->position);
    append(matePositions, firstGoodWin->matePosition);
    add(avgPerSampleSupport, firstGoodWin->perSampleSupport);
    double avgTotalSupport = firstGoodWin->totalSupport;
    unsigned winCount = 1;
    unsigned callCount = 1;

    Iterator<String<JunctionCall>, Standard >::Type it = firstGoodWin + 1;
    while (true)
    {
        if (similar(*currentAnchor, *it, stddev))
        {
            avgTotalSupport += it->totalSupport;
            add(avgPerSampleSupport, it->perSampleSupport);
            //std::cout << it << " " << it->windowPosition << " " << it->position << " " << it->matePosition << " " << it->perSampleSupport[0] << std::endl;
            appendValue(positions, it->position);
            appendValue(matePositions, it->matePosition);
            ++winCount;
            invalidateJunctionCall(*it);
            if (it == last)
            {
                if (!empty(positions))
                    mergeJunctionCallsRange(currentAnchor,
                                            it,
                                            positions,
                                            matePositions,
                                            avgPerSampleSupport,
                                            avgTotalSupport,
                                            genotypes,
                                            callCount,
                                            winCount);
                break;
            }
        }
        else
        {
            if (winCount != 1 && !empty(positions))
                mergeJunctionCallsRange(currentAnchor,
                                        it,
                                        positions,
                                        matePositions,
                                        avgPerSampleSupport,
                                        avgTotalSupport,
                                        genotypes,
                                        callCount,
                                        winCount);
            else
                invalidateJunctionCall(*currentAnchor);

            currentAnchor = it;
        }
        if (it != last)
        {
            ++it;
        }
        else
        {
            if (winCount == 1)
            {
                invalidateJunctionCall(*currentAnchor);
            }
            break;
        }
    }
    reserve(tmp, callCount, Exact());
    it = firstGoodWin;
    for (; it <= last; ++it)
    {
        if (isValid(*it))
            appendValue(tmp, *it);
    }
    move(junctionCalls, tmp);
    return true;
}
// Overload for working with translocation junctions
inline bool mergeJunctionWindows(String<JunctionCall> & junctionCalls, const double & stddev, TranslocJunctions t)
{
    (void) t;
    if (length(junctionCalls) < 1u)
    {
        return false;
    }
    else if (length(junctionCalls) == 1u)
    {
       refineMergedTranslocCall(*begin(junctionCalls));
       return true;
    }
    std::cout << "\nSorting " << length(junctionCalls) << " translocCalls" << std::endl;
    std::cout << "Before sorting:"<< std::endl;
    for (auto it = begin(junctionCalls); it != end(junctionCalls); ++it)
    {
        it->print();
        std::cout << "\nFR-Pos:\t";
        for (auto itt = begin(it->perSampleFRPositions); itt != end(it->perSampleFRPositions); ++itt)
        {
            std::cout << *itt << ", ";
        }
        std::cout << std::endl;
        std::cout << "RF-Pos:\t";
        for (auto itt = begin(it->perSampleRFPositions); itt != end(it->perSampleRFPositions); ++itt)
        {
            std::cout << *itt << ", ";
        }
        std::cout << std::endl;
        std::cout << "SupportCount:\t";
        for (auto itt = begin(it->perSampleTranslocSupport); itt != end(it->perSampleTranslocSupport); ++itt)
        {
            std::cout << *itt << ", ";
        }
        std::cout << std::endl;
    }
    std::sort(begin(junctionCalls), end(junctionCalls), translocOrder);
    std::cout << "\nAfter sorting:"<< std::endl;
    for (auto it = begin(junctionCalls); it != end(junctionCalls); ++it)
    {
        it->print();
        std::cout << std::endl;
    }
    std::cout << std::endl;

    Iterator<String<JunctionCall> >::Type currentAnchor = begin(junctionCalls, Standard());
    Iterator<String<JunctionCall> >::Type nextAnchorCandidate = currentAnchor;
    unsigned callCount = 1;

    for (Iterator<String<JunctionCall>, Standard >::Type it = currentAnchor + 1 ;it != end(junctionCalls); ++it)
    {
        if (isValid(*it))   // Skip calls that have already been merged
        {
            if (similar(*currentAnchor, *it, stddev, TranslocJunctions()))
            {
                pairwiseTranslocMerge(currentAnchor, it); // Merge and invalidate it
            }
            else if (currentAnchor->mateChromosome == it->mateChromosome &&  // can there still be merge candidates?
                     currentAnchor->position + 2 * stddev >= it->position)
            {
                if (nextAnchorCandidate == currentAnchor)   // Set new anchor candiate if not yet set.
                {
                    nextAnchorCandidate = it;
                }
            }
            else    // No more possible merge candidates for this anchor. Move it.
            {
                ++callCount;
                if (nextAnchorCandidate > currentAnchor)
                {
                    currentAnchor = nextAnchorCandidate;
                    it = currentAnchor; // Move iterator back to the new anchor to not miss any merge candidates.
                }
                else
                {
                    currentAnchor = it;
                }
            }
        }
    }
    String<JunctionCall> tmp;
    reserve(tmp, callCount, Exact());
    for (Iterator<String<JunctionCall>, Standard >::Type it = begin(junctionCalls); it != end(junctionCalls); ++it)
    {
        if (isValid(*it))
        {
            refineMergedTranslocCall(*it);
            appendValue(tmp, *it);
        }
    }
    move(junctionCalls, tmp);

    std::cout << "\nAfter merging:"<< std::endl;
    for (auto it = begin(junctionCalls); it != end(junctionCalls); ++it)
    {
        it->print();
                std::cout << "\nFR-Pos:\t";
        for (auto itt = begin(it->perSampleFRPositions); itt != end(it->perSampleFRPositions); ++itt)
        {
            std::cout << *itt << ", ";
        }
        std::cout << std::endl;
        std::cout << "RF-Pos:\t";
        for (auto itt = begin(it->perSampleRFPositions); itt != end(it->perSampleRFPositions); ++itt)
        {
            std::cout << *itt << ", ";
        }
        std::cout << std::endl;
        std::cout << "SupportCount:\t";
        for (auto itt = begin(it->perSampleTranslocSupport); itt != end(it->perSampleTranslocSupport); ++itt)
        {
            std::cout << *itt << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return true;
}
// Overload for application on tuple of differently oriented junctionWindows.
inline bool mergeJunctionWindows(Tuple<String<JunctionCall>, 4> & junctionCalls, const double & stddev)
{
    bool rf = mergeJunctionWindows(junctionCalls[0], stddev);
    bool ff = mergeJunctionWindows(junctionCalls[1], stddev);
    bool rr = mergeJunctionWindows(junctionCalls[2], stddev);
    bool tt = false; // mergeJunctionWindows(junctionCalls[3], stddev, TranslocJunctions());
    return (rf || ff || rr || tt);
}
#endif /* DETECT_JUNCTION_POPDEL_CALL_H_ */
