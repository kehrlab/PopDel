#ifndef INCLUDE_SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_UTILS_H_
#define INCLUDE_SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_UTILS_H_

#include <seqan/bam_io/bam_alignment_record.h>
#include <seqan/hts_io/hts_file.h>
#include <unordered_set>

namespace seqan {


CharString inline
getContigName(int32_t const rID, HtsFile const & file)
{
    SEQAN_ASSERT(file.hdr);

    if (rID < 0 || rID >= sam_hdr_nref(file.hdr))
        return CharString("");

    return CharString(file.hdr->target_name[rID]);
}
std::string inline
getContigNameAsStdString(int32_t const rID, HtsFile const & file)
{
    SEQAN_ASSERT(file.hdr);

    if (rID < 0 || rID >= sam_hdr_nref(file.hdr))
        return std::string("");

    return std::string(file.hdr->target_name[rID]);
}
// =======================================================================================
// Function getContigNames()
// =======================================================================================
// Get a String of all contig names
void inline
getContigNames(String<String<char> > & contigNames, HtsFile const & file)
{
    const int32_t nContigs = sam_hdr_nref(file.hdr);
    if (nContigs <= 0)
        return;

    reserve(contigNames, nContigs, Exact());
    for (int32_t rID = 0; rID < nContigs; ++rID)
        appendValue(contigNames, getContigName(rID, file));
}
void inline
getContigNames(std::unordered_set<std::string> & contigNames, HtsFile const & file)
{
    const int32_t nContigs = sam_hdr_nref(file.hdr);
    if (nContigs <= 0)
        return;

    contigNames.reserve(nContigs);
    for (int32_t rID = 0; rID < nContigs; ++rID)
        contigNames.emplace(getContigNameAsStdString(rID, file));
}
// =======================================================================================
// Function getContigNameToIDMap()
// =======================================================================================
// Create a map of contig names to rID.
void inline
getContigNameToIDMap(std::map<CharString, int32_t > & contigNames, HtsFile const & file)
{
    const int32_t nContigs = sam_hdr_nref(file.hdr);
    if (nContigs <= 0)
        return;

    for (int32_t rID = 0; rID < nContigs; ++rID)
        contigNames.emplace(std::make_pair(getContigName(rID, file), rID));
}
// =======================================================================================
// Function getContigLength()
// =======================================================================================
// Return length of the contig with the given rID
int inline
getContigLength(int32_t const rID, HtsFile const & file)
{
    SEQAN_ASSERT(file.hdr);

    if (rID < 0 || rID >= sam_hdr_nref(file.hdr))
        return -1;

    return file.hdr->target_len[rID];
}
// =======================================================================================
// Function getContigLengths()
// =======================================================================================
// Get lengths of the all contigs.
void inline
getContigLengths(String<int32_t> & contigLengths, HtsFile const & file)
{
    //const int32_t nContigs = file.hdr->n_targets;
    const int32_t nContigs = sam_hdr_nref(file.hdr);
    reserve(contigLengths, nContigs, Exact());
    for (int32_t rID = 0; rID < nContigs; ++rID)
        appendValue(contigLengths, getContigLength(rID, file));
}

} // namespace seqan

#endif // INCLUDE_SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_UTILS_H_
