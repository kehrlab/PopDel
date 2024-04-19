#ifndef PROFILE_PARAMETER_PARSING_POPDEL_H_
#define PROFILE_PARAMETER_PARSING_POPDEL_H_

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "bam_qual_req_popdel.h"
#include "../insert_histogram_popdel.h"
#include "../utils_popdel.h"

using namespace seqan;

// -----------------------------------------------------------------------------
// struct PopDelProfileParameters
// -----------------------------------------------------------------------------
// Struct holding all program parameters for computation of the insert size profile.
struct PopDelProfileParameters
{
    typedef std::map<CharString, unsigned> TReadGroups;
    CharString bamfile;                                   // Input file
    CharString outfile;                                   // Output file
    TReadGroups readGroups;                               // Map of read group ID and their number of appearance
    CharString intervalFile;                              // File of regular genomic intervals for parameter estimation.
    String<Histogram> histograms;                         // Insert size distributions per read group.
    String<unsigned> sampleSize;                          // Counts the number of sampled reads per read group.
    unsigned maxDeletionSize;                             // Iteration parameters. Maximum size of a deletion.
    unsigned windowSize;                                  // Size of the window.
    unsigned windowShift;                                 // Bases by which the window is shifted each iteration.
    unsigned minSampling;                                 // Min number of read pairs for hists per read group.
    String<GenomicRegion> rois;                           // Genomic intervals to iterate with the window iterator.
    BamQualReq qualReq;                                   // Quality requirements of the reads.
    unsigned indexRegionSize;                             // Size of regions for which file offsets are stored in profile index.
    bool mergeRG;
    CharString referenceVersion;
    CharString referenceFile;
    bool uncompressed;

    PopDelProfileParameters() :
    outfile("*BAM/CRAM-FILE*.profile"),
    intervalFile(""),
    maxDeletionSize(10000),
    windowSize(256),
    windowShift(0),
    minSampling(50000),
    indexRegionSize(10000),
    mergeRG(false),
    referenceVersion("GRCh38"),
    referenceFile(""),
    uncompressed(false)
    {}
};
// ---------------------------------------------------------------------------------------
// Function setHiddenOptions()
// ---------------------------------------------------------------------------------------
// Reveal/hide certain advanced options.
void setHiddenOptions(ArgumentParser & parser, bool hide, const PopDelProfileParameters &)
{
    hideOption(parser, "ir", hide);
    hideOption(parser, "f",  hide);
    hideOption(parser, "F",  hide);
    hideOption(parser, "mq", hide);
    hideOption(parser, "u",  hide);
    hideOption(parser, "R",  hide);
    hideOption(parser, "s",  hide);
    hideOption(parser, "x",  hide);
}
// ---------------------------------------------------------------------------------------
// Function addHiddenOptions()
// ---------------------------------------------------------------------------------------
// Add advanced options to the parser, which are only visible in the full help.
void addHiddenOptions(ArgumentParser & parser, const PopDelProfileParameters & params)
{
    addOption(parser, ArgParseOption("f",  "flags-set",        "Only use reads with all bits of NUM set in the bam flag.",                   ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("F",  "flags-unset",      "Only use reads with all bits of NUM unset in the bam flag.",                 ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("ir", "index-region-size", "Size of the index region intervals.",                                       ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("mq", "min-mapping-qual", "Only use reads with a mapping quality above NUM.",                           ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("R", "reference-file",    "FASTA file of the reference genome. Only required for CRAM files whose header entries point to the wrong file.", ArgParseArgument::STRING, "FILE"));
    addOption(parser, ArgParseOption("s",  "min-align-score",  "Only use reads with an alignment score relative to read length above NUM.",  ArgParseArgument::DOUBLE,  "NUM"));
    addOption(parser, ArgParseOption("u",  "min-unclipped",    "Only use reads of which at least NUM bases are not clipped.",                ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("x",  "uncompressed",    "Don't use gzip compression for profile."));

    setDefaultValue(parser, "f",  params.qualReq.flagsSet);
    setDefaultValue(parser, "F",  params.qualReq.flagsUnset);
    setDefaultValue(parser, "mq", params.qualReq.minMappingQual);
    setDefaultValue(parser, "u",  params.qualReq.minUnclippedLength);
    setDefaultValue(parser, "s",  params.qualReq.minRelAlignScore);
    setDefaultValue(parser, "ir",  params.indexRegionSize);
    setHiddenOptions(parser, true, params);
}
// -----------------------------------------------------------------------------
// Function setupParser()
// -----------------------------------------------------------------------------
// Set up the argument parser.
void setupParser(ArgumentParser & parser, const PopDelProfileParameters & params)
{
    setShortDescription(parser, "Profile creation from BAM/CRAM file");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM/CRAM-FILE\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM/CRAM-FILE\\fP \\fICHROM:BEGIN-END\\fP [\\fICHROM:BEGIN-END\\fP ...]");
    addDescription(parser, "Iterates over (user definied regions of) a BAM or CRAM file in tiling windows of 256 bp. For "
    " each window, PopDel promotes all read pairs whose ends pass the quality checks."
    " PopDel saves their insert size deviation form the mean together with their position in \'*BAM/CRAM-FILE*.profile\',"
    " together with the insert sizes distribution of each read group."
    " Only insert sizes up to a maximum length (option --max-deletion-size) are considered.");
    // Require a bam/cram file as argument and optionally a list of genomic intervals.
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "INPUT", true));
    // Add visible options.
    addOption(parser, ArgParseOption("H", "fullHelp",          "Displays full list of options."));
    addSection(parser, "PopDel profile options");
    addOption(parser, ArgParseOption("r", "reference",         "Reference genome version used for the mapping. Not used when using custom sampling intervals (option '-i'). One of 'T2T', 'GRCh38', 'GRCh37', 'hg38', 'hg19' (case-insensitive).", ArgParseArgument::STRING, "REF"));
    addOption(parser, ArgParseOption("d", "max-deletion-size", "Maximum size of deletions.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("i", "intervals",         "File with genomic intervals for parameter estimation instead of default intervals (see README). One closed interval per line, formatted as \'CHROM:START-END\', 1-based coordinates.", ArgParseArgument::INPUT_FILE, "FILE"));
    addOption(parser, ArgParseOption("mrg","merge-read-groups","Merge all read groups of the sample. Only advised if they share the same properties!"));
    addOption(parser, ArgParseOption("n", "min-read-num",      "Minimum number of read pairs for parameter estimation (per read group)", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("o", "out",               "Output file name.", ArgParseArgument::OUTPUT_FILE, "FILE"));
    setDefaultValue(parser, "d", params.maxDeletionSize);
    setDefaultValue(parser, "n", params.minSampling);
    setDefaultValue(parser, "o", params.outfile);
    setDefaultValue(parser, "r", "GRCh38");
    addHiddenOptions(parser, params);                            // Add options that are visible only in the full help.
}
// ---------------------------------------------------------------------------------------
// Function getArgumentValues()
// ---------------------------------------------------------------------------------------
// Get the path to input file and parse the genomic intervals (if given via command line)
// Calls parseIntervals and mergeOverlappingIntervals.
void getArgumentValues(PopDelProfileParameters & params, ArgumentParser & parser)
{
    std::vector<std::string> arguments = getArgumentValues(parser, 0);
    params.bamfile = arguments[0];
    if (length(arguments) > 1)
    {
        arguments.erase(arguments.begin());
        parseIntervals(params.rois, arguments, params.bamfile);
//        expandIntervals(params.intervals, params.maxDeletionSize);
        mergeOverlappingIntervals(params.rois, params.maxDeletionSize);    //Merge ROIs from command line.
    }
}
// ---------------------------------------------------------------------------------------
// Function parseReferenceGenome()
// ---------------------------------------------------------------------------------------
// Parse the version of the reference genome. Ignore case.
inline void parseReferenceGenome(CharString & refVersion, const ArgumentParser & parser)
{
    String<char> ref;
    if(!getOptionValue(ref, parser, "reference"))
        return;     // Option not set -> Assume GRCh38.

    if (isSet(parser, "intervals"))
    {
        printStatus("WARNING: Custom sampling intervals have been defined (option '-i') and a reference"
                    " genome has been specified (option '-r'). Option '-r' will be ignored.");
        return;
    }
    refVersion = ref;
    toLower(refVersion);
    if (refVersion != "grch38" && refVersion != "hg38" && refVersion != "grch37" && refVersion != "hg19" && refVersion != "t2t")
    {
        std::ostringstream msg;
        msg << "[PopDel] If defined, reference genome version (-r, --reference) must be one of 'T2T', 'GRCh38', 'GRCh37',"
               " 'hg38', 'hg19' (case-insensitive). Got '" << ref << "' instead. Terminating.";
        SEQAN_THROW(ParseError(toCString(msg.str())));
    }
}
// ---------------------------------------------------------------------------------------
// Function getParameterValues()
// ---------------------------------------------------------------------------------------
// Get all parameter values from the parser.
void getParameterValues(PopDelProfileParameters & params, const ArgumentParser & parser)
{
    getOptionValue(params.maxDeletionSize,            parser, "max-deletion-size");
    getOptionValue(params.minSampling,                parser, "min-read-num");
    getOptionValue(params.outfile,                    parser, "out");
    getOptionValue(params.intervalFile,               parser, "intervals");
    getOptionValue(params.qualReq.flagsSet,           parser, "flags-set");
    getOptionValue(params.qualReq.flagsUnset,         parser, "flags-unset");
    getOptionValue(params.qualReq.minMappingQual,     parser, "min-mapping-qual");
    getOptionValue(params.qualReq.minUnclippedLength, parser, "min-unclipped");
    getOptionValue(params.qualReq.minRelAlignScore,   parser, "min-align-score");
    getOptionValue(params.indexRegionSize,            parser, "index-region-size");
    getOptionValue(params.referenceFile,              parser, "reference-file");
    parseReferenceGenome(params.referenceVersion,     parser);
    params.mergeRG = isSet(                           parser, "merge-read-groups");
    params.uncompressed = isSet(                      parser, "uncompressed");
}

#endif /* PROFILE_PARAMETER_PARSING_POPDEL_H_ */
