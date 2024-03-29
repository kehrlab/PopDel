#ifndef VCFOUT_POPDEL_H_
#define VCFOUT_POPDEL_H_

#include <chrono>

#include <seqan/vcf_io.h>

#include "../utils_popdel.h"
#include "parameter_parsing_popdel_call.h"

using namespace seqan;

VcfHeader buildVcfHeader(const String<CharString> & contigNames, const String<int32_t> & contigLengths)
{
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::ostringstream date;
    char timeBuffer[24];
    if (0 < std::strftime(timeBuffer, sizeof(timeBuffer), "[%Y%m%d] ", std::localtime(&now_c)))
        date << timeBuffer;
    String<char> source = "PopDel-V";
    append(source, VERSION);
    VcfHeader header;
    appendValue(header, VcfHeaderRecord("fileformat", "VCFv4.3"));
    appendValue(header, VcfHeaderRecord("fileDate", date.str()));
    appendValue(header, VcfHeaderRecord("source", source));
    for (unsigned i = 0; i < length(contigNames); ++i)
    {
        CharString s;
        append(s, "<ID=");
        append(s, contigNames[i]);
        append(s, ",length=");
        append(s, std::to_string(contigLengths[i]));
        append(s, ">");
        appendValue(header, VcfHeaderRecord("contig", s));
    }
    appendValue(header, VcfHeaderRecord("INFO", "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=SVMETHOD,Number=1,Type=String,Description=\"Approach used to detect the structural variant\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=LR,Number=1,Type=String,Description=\"Log-Likelihood ratio that the test is correct\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=YIELD,Number=1,Type=Float,Description=\"Fraction of genotyped samples\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=SWIN,Number=1,Type=Integer,Description=\"Number of significant windows merged into this variant\">"));
    appendValue(header, VcfHeaderRecord("FILTER", "<ID=LowLR,Description=\"Likelihood ratio below threshold\">"));
    //appendValue(header, VcfHeaderRecord("FILTER", "<ID=highCov,Description=\"High coverage region\">")); // TODO. Add Filter
    appendValue(header, VcfHeaderRecord("FILTER", "<ID=missingSamples,Description=\"Too many samples not genotyped\">"));
    appendValue(header, VcfHeaderRecord("FILTER", "<ID=allRefGT,Description=\"All samples genotyped as homozygous reference\">"));
    appendValue(header, VcfHeaderRecord("FILTER", "<ID=CSWin,Description=\"Low fraction of significant windows\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality. Difference of the best and second-best PL\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=LAD,Number=3,Type=Integer,Description=\"Likelihood derived allelic depth: Count of read-pairs supporting REF, ambiguous, ALT\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=DAD,Number=5,Type=Integer,Description=\"Distribution derived allelic depth: Count of read-pairs supporting REF only, REF and ALT, neither(between the histograms), ALT only, bigger than ALT\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=FL,Number=2,Type=Integer,Description=\"Window of the first and last read active in this window\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=FLD,Number=1,Type=Integer,Description=\"Distance between first and last window\">"));
    return header;
}
// Calculate the GQ from the PL values.
inline unsigned calculateGQ(const Triple<unsigned> & pl)
{
    if (sum(pl) == 0u)
        return 0;

    if (pl.i1 == 0u)
        return std::min(pl.i2, pl.i3);
    else if (pl.i2 == 0u)
        return std::min(pl.i1, pl.i3);
    else
        return std::min(pl.i1, pl.i2);
}
// Create the format string.
inline void buildFormatString(VcfRecord & record,
                              const Triple<unsigned> & phredLikelihoods,
                              const Triple<unsigned> & lads,
                              const Dad & dads,
                              const Pair<unsigned> & mappDist)
{
    unsigned gq = std::min(calculateGQ(phredLikelihoods), 255u);
    std::ostringstream formatString;
    if (gq == 0u)
    {
        formatString << "./.:0,0,0";
    }
    else
    {
        if (phredLikelihoods.i3 == 0)                                       // Hom. deletion.
            formatString << "1";
        else
            formatString << "0";
        formatString << "/";
        if (phredLikelihoods.i1 != 0)
            formatString << "1";                                            // Het. or hom. deletion.
        else
            formatString << "0";                                            // No deletion.
        formatString << ":"
                     << std::min(phredLikelihoods.i1, 255u) << ",";         // use 255 as upper bound for phred-score.
        formatString << std::min(phredLikelihoods.i2, 255u) << ",";
        formatString << std::min(phredLikelihoods.i3, 255u);
    }
    formatString << ":" << gq;
    formatString << ":" << lads.i1 << "," << lads.i2 << "," << lads.i3;
    formatString << ":" << dads.ref << "," << dads.both << "," << dads.between << "," << dads.alt << "," << dads.right;
    formatString << ":" << mappDist.i1 << "," << mappDist.i2;
    formatString << ":" << mappDist.i2 - mappDist.i1;
    appendValue(record.genotypeInfos, formatString.str());
}
inline double getYield(const Call & call)
{
    unsigned n = length(call.gtLikelihoods);
    double absYield = n;
    for (Iterator<const String<Triple<unsigned> > >::Type it = begin(call.gtLikelihoods); it != end(call.gtLikelihoods); ++it)
        if (sum(*it) == 0u)
            --absYield;
    return absYield / n;
}
inline String<char> setFilterString(const Call & call)
{
    if (checkAllPass(call))
    {
        return "PASS";
    }
    else
    {
        bool oneSet = false;
        String<char> filter("");
        if (!checkLRPass(call))
        {
            filter = "lowLR";
            oneSet = true;
        }
        if (!checkCoveragePass(call))
        {
            if (oneSet)
            {
                appendValue(filter, ';');
            }
            append(filter, "highCov");
            oneSet = true;
        }
        if (!checkSamplePass(call))
        {
            if (oneSet)
            {
                appendValue(filter, ';');
            }
            append(filter, "missingSamples");
            oneSet = true;
        }
        if (!checkGT0Pass(call))
        {
            if (oneSet)
            {
                appendValue(filter, ';');
            }
            append(filter, "allRefGT");
            oneSet = true;
        }
        if (!checkRelWinCovPass(call))
        {
            if (oneSet)
            {
                appendValue(filter, ';');
            }
            append(filter, "CSWin");
            //oneSet = true;  // Uncomment if more filters are added.
        }
        return filter;
    }
}
inline unsigned getEnd(const unsigned pos, const unsigned len)
{
    return pos + len;
}
inline VcfRecord buildRecord(const QuantileMap& quantileMap, const Call& call, int32_t rID)
{
    VcfRecord record;
    record.rID = rID;
    if (call.position > 1)
        record.beginPos = call.position - 1;
    else
        record.beginPos = call.position;
    record.id = ".";
    record.ref = "N";
    record.alt = "<DEL>";
    double errorProb = 1 - quantileMap.getProbabilty(call.likelihoodRatio);
    record.qual = errorProb == 0 ? 100 : _round(-10 * log10(errorProb));
    record.filter = setFilterString(call);
    std::ostringstream info;
    info << "IMPRECISE;"
         << "SVLEN=" << - static_cast<int>(call.deletionLength) << ";"
         << "END="<< getEnd(call.position, call.deletionLength) << ";"
         << "SVTYPE=DEL;"
         << "AF=" << call.frequency << ";"
         << "LR=" << call.likelihoodRatio << ";"
         << "SVMETHOD=PopDelv" << VERSION << ";"
         << "YIELD=" << getYield(call) << ";"
         << "SWIN=" << call.significantWindows;
    record.info =  info.str();
    record.format = "GT:PL:GQ:LAD:DAD:FL:FLD";
    for (unsigned i = 0; i < length(call.gtLikelihoods); ++i)
        buildFormatString(record, call.gtLikelihoods[i], call.lads[i], call.dads[i], call.firstLast[i]);
    return record;
}
// Hold everything necessary for the VCF-output.
struct VcfOutputBundle
{
    std::fstream out;
    VcfFileOut vcfOut;
    unsigned rID;                                                  // Keeps track of the ID for the current chromosome.
    bool contigNameLock;                                           // True if the current contig already has been added.

    VcfOutputBundle(const CharString & outfile,
                    const String<CharString> & inputFiles,
                    const String<CharString> & contigNames,
                    const String<int32_t> & contigLengths)
    {
        out.open(toCString(outfile), std::ios::out);
        if (!out.good())
            SEQAN_THROW(FileOpenError(toCString(outfile)));
        open(vcfOut, out, Vcf());
        appendSampleNames(inputFiles);
        VcfHeader vcfHeader = buildVcfHeader(contigNames, contigLengths);
        writeHeader(vcfOut, vcfHeader);
        out << std::flush;
        rID = 0;
        contigNameLock = false;
    }

    inline void flush()
    {
        out << std::flush;
    }

    inline void lockContigName()
    {
        SEQAN_ASSERT(!contigNameLock);
        contigNameLock = true;
    }

    inline void unlockContigName()
    {
        if (contigNameLock)
        {
            ++rID;
            contigNameLock = false;
        }
    }
    // Only works if contigNameLock == false !
    // Call unlockContigName() to remove lock.
    inline void appendContigName(const CharString & name)
    {
        if (!contigNameLock)
        {
            appendValue(contigNames(context(vcfOut)), name);
            lockContigName();
        }
    }

    inline void appendSampleNames(const String<CharString> & inputFiles)
    {
        for (unsigned i = 0; i < length(inputFiles); ++i)
        {
            appendValue(sampleNames(context(vcfOut)), getSampleName(inputFiles[i]));
        }
    }

    inline void outputRecord(VcfRecord & record)
    {
        writeRecord(vcfOut, record);
    }
};
// Iterate over all re-genotyped calls and create their VCF-records.
inline void writeRegenotypedCalls(VcfOutputBundle & vcfOutput,
                                  const String<Call> & reGenotypedCalls,
                                  const PopDelCallParameters & params)
{
    typedef Iterator<const String<Call>, Rooted>::Type TBufferIt;
    for (TBufferIt it = begin(reGenotypedCalls); it != end(reGenotypedCalls); ++it)
    {
        if (isValidCall((*it)))
        {
            if (params.outputFailed || checkAllPass(*it))
            {
                VcfRecord record = buildRecord(params.quantileMap, *it, vcfOutput.rID);
                vcfOutput.outputRecord(record);
            }
        }
    }
}
#endif /* VCFCOUT_POPDEL_H_ */
