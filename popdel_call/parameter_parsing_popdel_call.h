#ifndef PARAMETER_PARSING_POPDEL_CALL_H_
#define PARAMETER_PARSING_POPDEL_CALL_H_

#include <seqan/arg_parse.h>

#include "../insert_histogram_popdel.h"
#include "utils_popdel.h"

using namespace seqan;
// -----------------------------------------------------------------------------
// struct QuantileMap
// -----------------------------------------------------------------------------
//Maps a log likelohood ratio to the probability that the call is correct (the logLR has to be scaled correctly)
struct QuantileMap
{
    std::map<double, double> m;
    // Key: LogLR, value: quantile.
    // unscaled keys calculated with qchisq(value, df=1) in Rstudio and scaled
    // Values selected as probabilities of phred scores 0 to 50
    // Code for generation in R:
    // options(digits=15)
    // for (i in 0:100)
    //     cat("m.emplace(", qchisq(1-10^(-i/10), df=1), " / 2 - p, ", 1-10^(-i/10), ");\n", sep = '')
    QuantileMap(double prior)
    {
        double p = log(prior / (1 - prior));
        m.emplace(0, 0);
        m.emplace(0.0679615365255401 / 2 -p,0.205671765275719);
        m.emplace(0.230764766661193 / 2 - p,0.369042655519807);
        m.emplace(0.452421556097698 / 2 - p,0.498812766372728);
        m.emplace(0.714036104152315 / 2 - p,0.601892829446503);
        m.emplace(1.00448471501902 / 2 - p, 0.683772233983162);
        m.emplace(1.31668063342863 / 2 - p, 0.748811356849042);
        m.emplace(1.64583886397469 / 2 - p, 0.800473768503112);
        m.emplace(1.98858107774449 / 2 - p, 0.841510680753889);
        m.emplace(2.34243609502603 / 2 - p, 0.874107458820583);
        m.emplace(2.70554345409542 / 2 - p, 0.9);
        m.emplace(3.07646857913928 / 2 - p, 0.920567176527572);
        m.emplace(3.45408274293754 / 2 - p, 0.936904265551981);
        m.emplace(3.83748240184564 / 2 - p, 0.949881276637273);     // ~95%-quantile
        m.emplace(4.22593339019369 / 2 - p, 0.96018928294465);
        m.emplace(4.61883133337202 / 2 - p, 0.968377223398316);
        m.emplace(5.01567294625387 / 2 - p, 0.974881135684904);
        m.emplace(5.41603482049954 / 2 - p, 0.980047376850311);
        m.emplace(5.81955747770787 / 2 - p, 0.984151068075389);
        m.emplace(6.22593319775977 / 2 - p, 0.987410745882058);
        m.emplace(6.63489660102121 / 2 - p, 0.99);
        m.emplace(7.04621727098416 / 2 - p, 0.992056717652757);
        m.emplace(7.45969391025114 / 2 - p, 0.993690426555198);
        m.emplace(7.87514966368955 / 2 - p, 0.994988127663727);
        m.emplace(8.29242834051146 / 2 - p, 0.996018928294465);
        m.emplace(8.71139133617301 / 2 - p, 0.996837722339832);
        m.emplace(9.13191510451069 / 2 - p, 0.99748811356849);
        m.emplace(9.55388906647803 / 2 - p, 0.998004737685031);
        m.emplace(9.97721386826398 / 2 - p, 0.998415106807539);
        m.emplace(10.4017999212068 / 2 - p, 0.998741074588206);
        m.emplace(10.8275661706627 / 2 - p, 0.999);
        m.emplace(11.2544390521783 / 2 - p, 0.999205671765276);
        m.emplace(11.6823516018745 / 2 - p, 0.99936904265552);
        m.emplace(12.1112426945654 / 2 - p, 0.999498812766373);
        m.emplace(12.5410563882771 / 2 - p, 0.999601892829446);
        m.emplace(12.9717413578738 / 2 - p, 0.999683772233983);
        m.emplace(13.403250403671 / 2 - p,  0.999748811356849);
        m.emplace(13.8355400234795 / 2 - p, 0.999800473768503);
        m.emplace(14.2685700385079 / 2 - p, 0.999841510680754);
        m.emplace(14.7023032652319 / 2 - p, 0.999874107458821);
        m.emplace(15.1367052266236 / 2 - p, 0.9999);
        m.emplace(15.571743897219 / 2 - p,  0.999920567176528);
        m.emplace(16.007389477407 / 2 - p,  0.999936904265552);
        m.emplace(16.4436141930027 / 2 - p, 0.999949881276637);
        m.emplace(16.8803921167791 / 2 - p, 0.999960189282945);
        m.emplace(17.317699009203 / 2 - p,  0.999968377223398);
        m.emplace(17.7555121758597 / 2 - p, 0.999974881135685);
        m.emplace(18.193810339586 / 2 - p,  0.99998004737685);
        m.emplace(18.632573525518 / 2 - p,  0.999984151068075);
        m.emplace(19.0717829575225 / 2 - p, 0.999987410745882);
        m.emplace(19.5114209646663 / 2 - p, 0.99999);
        m.emplace(19.9514708965797 / 2 - p, 0.999992056717653);
        m.emplace(20.3919170468123 / 2 - p, 0.999993690426555);
        m.emplace(20.8327445832119 / 2 - p, 0.999994988127664);
        m.emplace(21.2739394844285 / 2 - p, 0.999996018928294);
        m.emplace(21.7154884822888 / 2 - p, 0.99999683772234);
        m.emplace(22.1573790090714 / 2 - p, 0.999997488113568);
        m.emplace(22.599599149316 / 2 - p,  0.999998004737685);
        m.emplace(23.0421375954437 / 2 - p, 0.999998415106808);
        m.emplace(23.4849836077724 / 2 - p, 0.999998741074588);
        m.emplace(23.9281269768795 / 2 - p, 0.999999);
        m.emplace(24.3715579896983 / 2 - p, 0.999999205671765);
        m.emplace(24.8152673978314 / 2 - p, 0.999999369042656);
        m.emplace(25.2592463886412 / 2 - p, 0.999999498812766);
        m.emplace(25.7034865580547 / 2 - p, 0.999999601892829);
        m.emplace(26.1479798870853 / 2 - p, 0.999999683772234);
        m.emplace(26.5927187165227 / 2 - p, 0.999999748811357);
        m.emplace(27.0376957272151 / 2 - p, 0.999999800473768);
        m.emplace(27.4829039228708 / 2 - p, 0.999999841510681);
        m.emplace(27.9283366047418 / 2 - p, 0.999999874107459);
        m.emplace(28.3739873627981 / 2 - p, 0.9999999);
        m.emplace(28.8198500523812 / 2 - p, 0.999999920567177);
        m.emplace(29.2659187856725 / 2 - p, 0.999999936904266);
        m.emplace(29.7121879123161 / 2 - p, 0.999999949881277);
        m.emplace(30.1586520172268 / 2 - p, 0.999999960189283);
        m.emplace(30.6053058897649 / 2 - p, 0.999999968377223);
        m.emplace(31.0521445329994 / 2 - p, 0.999999974881136);
        m.emplace(31.4991631422743 / 2 - p, 0.999999980047377);
        m.emplace(31.9463570799056 / 2 - p, 0.999999984151068);
        m.emplace(32.3937219226572 / 2 - p, 0.999999987410746);
        m.emplace(32.8412533514689 / 2 - p, 0.99999999);
        m.emplace(33.2889472854096 / 2 - p, 0.999999992056718);
        m.emplace(33.7367997370742 / 2 - p, 0.999999993690427);
        m.emplace(34.1848068388937 / 2 - p, 0.999999994988128);
        m.emplace(34.6329649593365 / 2 - p, 0.999999996018928);
        m.emplace(35.0812704794425 / 2 - p, 0.999999996837722);
        m.emplace(35.5297199531207 / 2 - p, 0.999999997488114);
        m.emplace(35.9783101801565 / 2 - p, 0.999999998004738);
        m.emplace(36.4270377410164 / 2 - p, 0.999999998415107);
        m.emplace(36.8758997410294 / 2 - p, 0.999999998741075);
        m.emplace(37.3248931065187 / 2 - p, 0.999999999);
        m.emplace(37.7740147322013 / 2 - p, 0.999999999205672);
        m.emplace(38.2232620681608 / 2 - p, 0.999999999369043);
        m.emplace(38.672632697292 / 2 - p,  0.999999999498813);
        m.emplace(39.1221230144588 / 2 - p, 0.999999999601893);
        m.emplace(39.5717319767953 / 2 - p, 0.999999999683772);
        m.emplace(40.0214554741213 / 2 - p, 0.999999999748811);
        m.emplace(40.4712914238083 / 2 - p, 0.999999999800474);
        m.emplace(40.9212389675878 / 2 - p, 0.999999999841511);
        m.emplace(41.371295018385 / 2 - p,  0.999999999874107);
        m.emplace(41.8214562029826 / 2 - p, 0.9999999999);
    }
    // Return the closest (lower) probability value stored in the map.
    // Return 1 if lr is bigger than the highest value in the map.
    inline double getProbabilty(double lr) const
    {
        std::map<double, double>::const_iterator it = m.upper_bound(lr);
        return (it == m.end() ? 1.0 : (--it)->second);
    }
};
// -----------------------------------------------------------------------------
// struct PopDelCallParameters
// -----------------------------------------------------------------------------
// Struct holding all parameters required by podel call.
struct PopDelCallParameters
{
    typedef std::map<CharString, unsigned> TReadGroups;
    String<CharString> inputFiles;                  // Path/filename of file containing the paths to the profiles.
    CharString histogramsFile;                      // Path/filename of input histogram file.
    CharString roiFile;                             // File containing the regions of interest.
    CharString maxLoadFile;                         // File containing the maximum load for each read group.
    CharString outfile;                             // Output path/filename.
    unsigned iterations;                            // Max number of iterations for determination of deletion length.
    double prior;                                   // Prior for a deletion.
    double minimumLikelihoodRatio;                  // Cutoff from Chi-squared Test.
    TReadGroups readGroups;                         // Map of read group ID and their order of appearance.
    TRGs    rgs;                                    // One string of indices for the read groups per sample.
    String<Histogram> histograms;                   // Insert size distributions per read group.
    bool modRgByFileName;                           // True if the file names should be addec to the RGIDs.
    bool representativeContigs;                     // True if only the first sample's contig names shall be used
    String<unsigned> minInitDelLengths;             // Minimum length required for a deletion at initialization.
    unsigned minLen;                                // Minimum length estimation for deletion during iteration.
    double minRelWinCover;                          // Minimum number for (#SignificantWindows * 30 / DelSize)
    unsigned windowSize;                            // Lenght of one window (in BP).
    unsigned windowShift;                           // Shift of the window per iterarion.   // TODO. Remove
    unsigned windowBuffer;                          // Number of windows to buffer.
    bool smoothing;                                 // Bool indicating if the histogramms shall be smoothed or not.
    std::vector<std::string> roiList;               // List of rois as the come directly from the parser.
    String<GenomicRegion> allRois;                  // Genomic intervals to iterate with the window iterator.
    Iterator<String<GenomicRegion>, Standard>::Type nextRoi; // Points to the next ROI in allRois.
    bool outputFailed;                              // Bool indicating that also FAILED calls shold be written to file.
    QuantileMap quantileMap;                        // Map for translating log LR into PHRED error-probabilities.
    unsigned fileCount;                             // Number of remaining files.
    unsigned sampleNum;                             // Total number of files/Samples.
    String<String<CharString> > contigNames;        // Names of the contigs. One String per file.
    String<String<int32_t> > contigLengths;         // Lenghts of above contigs. One String per file.
    unsigned maxDeletionSize;                       // Maximum size of a deletion.
    double minSampleFraction;                       // Minimum fraction of samples which has to have data in the window.
    double meanStddev;                              // The mean of all standard deviations of the histograms.
    bool windowWiseOutput;                          // Output each window after window-wise genotyping.
    String<unsigned> indexRegionSizes;
    unsigned defaultMaxLoad;                        // Default value for the maximum load.
    String<unsigned> maxLoad;                       // Maximum number of active reads per read group.
    unsigned pseudoCountFraction;                   // max(hist)/pseudoCountFraction = min value of hist
    bool uncompressedIn;                            // True if uncompressed profiles are read.
    bool somatic;                                   // If true, use genotype weights for cell populations instead of HWE

    PopDelCallParameters() :
    histogramsFile(""),
    outfile("popdel.vcf"),
    iterations(15),
    prior(0.0001),
    modRgByFileName(false),
    representativeContigs(true),
    minLen(maxValue<unsigned>()),
    minRelWinCover(0.5),
    windowSize(30),
    windowShift(0),
    windowBuffer(200000),
    smoothing(true),
    outputFailed(false),
    quantileMap(0.0001),
    fileCount(0),
    sampleNum(0),
    maxDeletionSize(10000),
    minSampleFraction(0.1),
    meanStddev(0.0),
    windowWiseOutput(false),
    defaultMaxLoad(100),
    pseudoCountFraction(500),
    uncompressedIn(false),
    somatic(false)
    {}
};
// ---------------------------------------------------------------------------------------
// Function setHiddenOptions()
// ---------------------------------------------------------------------------------------
// Reveal/hide certain advanced options.
void setHiddenOptions(ArgumentParser & parser, bool hide, const PopDelCallParameters &)
{
    hideOption(parser, "b", hide);
    hideOption(parser, "c", hide);
    hideOption(parser, "f", hide);
    hideOption(parser, "F", hide);
    hideOption(parser, "n", hide);
    hideOption(parser, "p", hide);
    hideOption(parser, "t", hide);
    hideOption(parser, "u", hide);
    hideOption(parser, "x", hide);
}
// ---------------------------------------------------------------------------------------
// Function addHiddenOptions()
// ---------------------------------------------------------------------------------------
// Add advanced options to the parser, which are only visible in the full help.
void addHiddenOptions(ArgumentParser & parser, const PopDelCallParameters & params)
{
    addOption(parser, ArgParseOption("b", "buffer-size",       "Number of buffered windows.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("c", "min-relative-window-cover", "Determines which fraction of a deletion has to be covered by significant windows.", ArgParseArgument::DOUBLE, "NUM"));
    addOption(parser, ArgParseOption("f", "pseudocount-fraction",   "The biggest likelihood of the background distribution will be divided by this value to determine the pseudocounts of the histogram. Bigger values boost the sensitivity for HET calls but also increase the chance of missclassifying HOMDEL or HOMREF as HET calls.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("F", "output-failed", "Also output calls which did not pass the filters."));
    addOption(parser, ArgParseOption("n", "no-regenotyping",   "Outputs every potential variant window without re-genotyping and merging."));
    addOption(parser, ArgParseOption("p", "prior-probability", "Prior probability of a deletion.",                  ArgParseArgument::DOUBLE, "NUM"));
    addOption(parser, ArgParseOption("t", "iterations",        "Number of iterations in EM for length estimation.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("u", "unsmoothed",        "Disable the smoothing of the insert size histogram."));
    addOption(parser, ArgParseOption("x", "uncompressed-in",   "Read uncompressed PopDel profiles"));

    setDefaultValue(parser, "b",  params.windowBuffer);
    setDefaultValue(parser, "f",  params.pseudoCountFraction);
    setDefaultValue(parser, "p",  params.prior);
    setDefaultValue(parser, "t", params.iterations);
    setDefaultValue(parser, "c", params.minRelWinCover);
    setMinValue(parser, "buffer-size", "10000");
    setMinValue(parser, "pseudocount-fraction", "50");
    setMinValue(parser, "prior-probability", "0.0");
    setMaxValue(parser, "prior-probability", "0.9999");
    setMinValue(parser, "min-relative-window-cover", "0.0");
    setMaxValue(parser, "min-relative-window-cover", "2.0");
    setHiddenOptions(parser, true, params);
}
// -----------------------------------------------------------------------------
// Function setupParser()
// -----------------------------------------------------------------------------
// Set up the argument parser.
void setupParser(ArgumentParser & parser, const PopDelCallParameters & params)
{
    setShortDescription(parser, "Population-wide deletion calling");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fPROFILE-LIST-FILE\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIPROFILE-FILE1\\fP \\fIPROFILE-FILE2\\fP [\\fIPROFILE-FILE3\\fP ...]");
    addDescription(parser, "Performs joint-calling of deletions using a list of profile-files"
                           " previously created using the \'popdel profile\' command."
                           " The input profiles are either specified directly as arguments or listed in PROFILE-LIST-FILE"
                           " (one filename per line).");
    // Require a single input file as argument or a list of profile files.
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "INPUT", true));
    // Add visible options.
    addOption(parser, ArgParseOption("H", "fullHelp",              "Displays full list of options."));
    addSection(parser, "PopDel call options");
    addOption(parser, ArgParseOption("A", "active-coverage-file",  "File with lines consisting of \"ReadGroup  maxCov\". If this value is reached no more new reads are loaded for this read group until the coverage drops again. Further, the sample will be excluded  from calling in high-coverage windows. A value of 0 disables the filter for the read group.", ArgParseArgument::INPUT_FILE, "FILE"));
    addOption(parser, ArgParseOption("a", "active-coverage",       "Maximum number of active read pairs (~coverage). This value is taken for all read groups that are not listed in \'active-coverage-file\'. Setting it to 0 disables the filter for all read groups that are not specified in \'active-coverage-file\'.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("d", "max-deletion-size",     "Maximum size of deletions.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("e", "per-sample-rgid",       "Internally modify each read group ID by adding the filename. This can be used if read groups across different samples have conflicting IDs."));
    addOption(parser, ArgParseOption("l", "min-init-length",       "Minimal deletion length at initialization of iteration. Default: \\fI4 * standard deviation\\fP.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("m", "min-length",            "Minimal deletion length during iteration. Default: \\fI95th percentile of min-init-lengths\\fP.", ArgParseArgument::INTEGER, "NUM"));
    addOption(parser, ArgParseOption("o", "out",                   "Output file name.", ArgParseArgument::OUTPUT_FILE, "FILE"));
    addOption(parser, ArgParseOption("r", "region-of-interest",    "Genomic region 'chr:start-end' (closed interval, 1-based index). Calling is limited to this region. "
                                                                   "Multiple regions can be defined by using the parameter -r multiple times.",ArgParseArgument::STRING, "REGION", "isList"));
    addOption(parser, ArgParseOption("R", "ROI-file",              "File listing one or more regions of interest, one region per line. "
                                                                   "See parameter -r.",                                                        ArgParseArgument::INPUT_FILE, "FILE"));
    addOption(parser, ArgParseOption("s",  "min-sample-fraction",  "Minimum fraction of samples which is required to have enough data in the window.", ArgParseArgument::DOUBLE, "NUM"));
    addOption(parser, ArgParseOption("C",  "cell-population",      "Use genotype weights for cell population instead of human population."));

    // Set default values for visible options.
    setDefaultValue(parser, "a", params.defaultMaxLoad);
    setDefaultValue(parser, "d", params.maxDeletionSize);
    setDefaultValue(parser, "o", params.outfile);
    setDefaultValue(parser, "s", params.minSampleFraction);
    // Set min and max values
    setMinValue(parser, "active-coverage", "0");
    setMinValue(parser, "min-sample-fraction", "0.0");
    setMaxValue(parser, "min-sample-fraction", "1.0");
    // Add options that are visible only in the full help.
    addHiddenOptions(parser, params);
}
// ---------------------------------------------------------------------------------------
// Function getArgumentValues()
// ---------------------------------------------------------------------------------------
// Get the path to input files.
void getArgumentValues(PopDelCallParameters & params, ArgumentParser & parser)
{
    params.inputFiles = getArgumentValues(parser, 0);
    if (length(params.inputFiles) == 1)
        loadFilenames(params.inputFiles);
    params.sampleNum = length(params.inputFiles);
    params.fileCount = params.sampleNum;
}
// ---------------------------------------------------------------------------------------
// Function getParameterValues()
// ---------------------------------------------------------------------------------------
// Get all parameter values from the parser.
void getParameterValues(PopDelCallParameters & params, ArgumentParser & parser)
{
    getOptionValue(params.outfile,              parser, "out");
    getOptionValue(params.prior,                parser, "prior-probability");
    getOptionValue(params.iterations,           parser, "iterations");
    getOptionValue(params.roiFile,              parser, "ROI-file");
    getOptionValue(params.maxLoadFile,          parser, "active-coverage-file");
    getOptionValue(params.minSampleFraction,    parser, "min-sample-fraction");
    getOptionValue(params.windowBuffer,         parser, "buffer-size");
    getOptionValue(params.minRelWinCover,       parser, "min-relative-window-cover");
    getOptionValue(params.defaultMaxLoad,       parser, "active-coverage");
    getOptionValue(params.pseudoCountFraction,  parser, "pseudocount-fraction");
    params.smoothing = !isSet(parser, "unsmoothed");
    if (isSet(parser, "min-init-length"))
    {
        resize(params.minInitDelLengths, 1);
        getOptionValue(params.minInitDelLengths[0], parser, "min-init-length");
    }
    if (isSet(parser, "no-regenotyping"))
    {
        params.windowWiseOutput = true;
    }
    if (isSet(parser, "min-length"))
    {
        getOptionValue(params.minLen, parser, "min-length");
    }
    if (isSet(parser, "region-of-interest"))
    {
        unsigned regionCount = getOptionValueCount(parser, "region-of-interest");
        resize(params.roiList, regionCount);
        for (unsigned i = 0; i < regionCount; ++i)
            getOptionValue(params.roiList[i], parser, "region-of-interest", i);
    }
    if (isSet(parser, "output-failed"))
    {
        params.outputFailed = true;
    }
    if (isSet(parser, "per-sample-rgid"))
    {
        params.modRgByFileName = true;
    }
    if (isSet(parser, "uncompressed-in"))
    {
        params.uncompressedIn = true;
    }
    if (isSet(parser, "cell-population"))
    {
        params.somatic = true;
    }
}

#endif /* PARAMETER_PARSING_POPDEL_CALL_H_ */
