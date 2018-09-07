#ifndef POPDEL_VIEW_PARAMETER_PARSING_POPDEL_H_
#define POPDEL_VIEW_PARAMETER_PARSING_POPDEL_H_

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>


// -----------------------------------------------------------------------------
// struct PopDelViewParameters
// -----------------------------------------------------------------------------
// Struct holding all program parameters for viewing an insert size profile.
struct PopDelViewParameters
{
    CharString infile;                                   // Input file

    bool writeHeader;
    bool writeOnlyHeader;
    bool writeHistograms;

    unsigned indexRegionSize;
    GenomicRegion region;

    PopDelViewParameters() :
        writeHeader(false), writeOnlyHeader(false), writeHistograms(false), indexRegionSize(10000)
    {}
};

// ---------------------------------------------------------------------------------------
// Function setHiddenOptions()
// ---------------------------------------------------------------------------------------
// Reveal/hide certain advanced options.
void setHiddenOptions(ArgumentParser & /*parser*/, bool /*hide*/, const PopDelViewParameters & /*params*/)
{
    //hideOption(parser, "i", hide);
}

// ---------------------------------------------------------------------------------------
// Function addHiddenOptions()
// ---------------------------------------------------------------------------------------
// Add advanced options to the parser, which are only visible in the full help.
void addHiddenOptions(ArgumentParser & parser, const PopDelViewParameters & params)
{
    //addOption(parser, ArgParseOption("f",  "flags-set",        "Only use reads with all bits of NUM set in the bam flag.",                                      ArgParseArgument::INTEGER, "NUM"));
    //setDefaultValue(parser, "f",  params.qualReq.flagsSet);
    setHiddenOptions(parser, true, params);
}

// -----------------------------------------------------------------------------
// Function setupParser()
// -----------------------------------------------------------------------------
// Set up the argument parser.
void setupParser(ArgumentParser & parser, const PopDelViewParameters & params)
{
    setShortDescription(parser, "viewing a popdel profile file");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    //addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIPROFILE-FILE\\fP");
    addUsageLine(parser, "\\fIPROFILE-FILE\\fP");
    addDescription(parser, "Displays a profile file in human readable format.");

    // Require a profile file as argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "INPUT"));

    // Add visible options.
    addOption(parser, ArgParseOption("H", "fullHelp", "Displays full list of options."));
    addOption(parser, ArgParseOption("e", "header", "Write the header."));
    addOption(parser, ArgParseOption("E", "onlyHeader", "Only write the header."));
    addOption(parser, ArgParseOption("i", "histograms", "Write insert size histograms."));
    addOption(parser, ArgParseOption("r", "region", "Limit view to this genomic region.", ArgParseArgument::STRING, "CHR:BEGIN-END"));

    // Add options that are visible only in the full help.
    addHiddenOptions(parser, params);
}

// ---------------------------------------------------------------------------------------
// Function getArgumentValues()
// ---------------------------------------------------------------------------------------
// Get the path to input file and parse the genomic intervals (if given via command line)
// Calls parseIntervals and mergeOverlappingIntervals.
void getArgumentValues(PopDelViewParameters & params, ArgumentParser & parser)
{
    getArgumentValue(params.infile, parser, 0);
}

// ---------------------------------------------------------------------------------------
// Function getParameterValues()
// ---------------------------------------------------------------------------------------
// Get all parameter values from the parser.
void getParameterValues(PopDelViewParameters & params, ArgumentParser & parser)
{
    getOptionValue(params.writeHeader, parser, "header");
    getOptionValue(params.writeOnlyHeader, parser, "onlyHeader");
    getOptionValue(params.writeHistograms, parser, "histograms");

    if (isSet(parser, "region"))
    {
        CharString region;
        getOptionValue(region, parser, "region");
        parse(params.region, region);
    }
}


#endif /* POPDEL_VIEW_PARAMETER_PARSING_POPDEL_H_*/
