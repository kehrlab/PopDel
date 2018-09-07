#ifndef BAM_QUAL_REQ_POPDEL_H_
#define BAM_QUAL_REQ_POPDEL_H_

#include <seqan/bam_io.h>

#include "../utils_popdel.h"

using namespace seqan;
// =======================================================================================
// Struct BamQualReq
// =======================================================================================
struct BamQualReq
{
    __uint32 flagsSet;                   // Required flags.
    __uint32 flagsUnset;                 // Forbidden falgs
    bool sameChrom;
    unsigned minMappingQual;             // minimum mapping quality of each sequence.
    unsigned minUnclippedLength;         // minimum length of sequences when unclipped.
    float minRelAlignScore;              // minimum alignment score relative to length of unclipped alignment.
    int alignmentScoreIndex;             // index of the alignment score in the BamTagsDict.

    BamQualReq() :
    flagsSet(33),           // read paired &  mate reverse strand
    flagsUnset(3868),       // read/mate unmapped, read reverse strand, not primary alignment, supplementary alignment
                            // read fails platform/vendor quality checks, read is PCR or optical duplicate
    sameChrom(true), minMappingQual(1), minUnclippedLength(50), minRelAlignScore(0.8), alignmentScoreIndex(0)
    {}
};

// ---------------------------------------------------------------------------------------
// Function getAlignScore()
// ---------------------------------------------------------------------------------------
// Return the extracted the alignment score             // TODO: Check if a method without BamTagsDict is faster.
inline unsigned getAlignScore(const CharString & tags, int & index)
{
    unsigned score = 0;
    BamTagsDict dict(tags);
    if (length(dict) <= static_cast<unsigned>(index))    // If the number of tags got too low, reset index.
    {
        index = 0;
        if(!findTagKey(index, dict, "AS"))
        {
            SEQAN_ASSERT_FAIL("Could not find Alignment score in record!");
            index = 0;
        }
        if (!extractTagValue(score, dict, index))
        {
            SEQAN_ASSERT_FAIL("Not a valid integer at pos %u!", index);
            index = 0;
        }
    }
    else
    {
        if (getTagKey(dict, index) != "AS")             // Check if Tag at index is still alignment score
        {
            index = 0;                                  // If not, search again and update index.
            if(!findTagKey(index, dict, "AS"))
            {
                SEQAN_ASSERT_FAIL("Could not find Alignment score in record!");
                index = 0;
            }
        }
        if (!extractTagValue(score, dict, index))
        {
            SEQAN_ASSERT_FAIL("Not a valid integer at pos %u!", index);
            index = 0;
        }
    }
    return score;
}
// =======================================================================================
// Function isReverseRead()
// =======================================================================================
// Performs checks if read is primary alignment and mate is not reverse -> Read is primary alignment and reverse read
inline bool isReverseRead(const BamAlignmentRecord & record)
{   //Originally: 2336 instead of 3884, but 3884 is more restrictive, so we now use that.
    return (record.flag & 3884) == 0;   // i.e if NOT(mate on reverse strand|| not prim. alignment || supp. alignment ||
                                        //            read/mate is unmapped || duplicate || read fails platform QC)
}
// =======================================================================================
// Function meetsRequirements()
// =======================================================================================
// Perform checks if all requirements are met by the record.
// Return true if all checks are passed, false otherwise.
inline bool meetsRequirements(const BamAlignmentRecord & record, BamQualReq & qualReq)
{
    if ((record.flag & qualReq.flagsUnset) != 0)                        // Check if the -F specified flags are unset
        return false;
    if ((record.flag & qualReq.flagsSet) != qualReq.flagsSet)           // Check if the -f specified flags are set.
        return false;
    if (qualReq.sameChrom && record.rID != record.rNextId)  //If required, check if the two reads align to the same chr.
        return false;
    if (record.mapQ < qualReq.minMappingQual)                           // Check mapping quality.
        return false;
    SEQAN_ASSERT_GT(length(record.cigar), 0u);
    unsigned fullLenght = length(record.seq);
    unsigned leftClip = 0;                                              // Number of clipped bases at 5' end.
    unsigned rightClip = 0;                                             // Number of clipped bases at 3' end
    if (record.cigar[0].operation == 'S')
        leftClip = record.cigar[0].count;
    if (record.cigar[length(record.cigar) - 1].operation == 'S')
        rightClip = record.cigar[length(record.cigar) - 1].count;
    unsigned unclippedLength = fullLenght - leftClip - rightClip;      // Number of bases left after clipping
    if (unclippedLength < qualReq.minUnclippedLength)
        return false;
    if (getAlignScore(record.tags, qualReq.alignmentScoreIndex) < unclippedLength * qualReq.minRelAlignScore)
        return false;                                                                      // Check the alignment score.
    return true;
}
// =======================================================================================
// Function readRecord()
// =======================================================================================
// Read a record from BAM-file and perform checks.
// Return true if all checks are passed, false otherwise.
inline bool readRecord(BamAlignmentRecord & record, BamFileIn & infile, BamQualReq & qualReq)
{
    while (!atEnd(infile))
    {
        try
        {
            readRecord(record, infile);
        }
        catch (Exception const & e)
        {
            std::ostringstream msg;
            msg << "WARNING: Could no read record in BAM-File. Skipping in file..." << std::endl;
            printStatus(msg);
            return false;
        }
        if (meetsRequirements(record, qualReq))
            return true;
    }
    return false;
}
#endif /* BAM_QUAL_REQ_POPDEL_H_ */
