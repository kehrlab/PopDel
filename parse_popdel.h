#ifndef PARSE_POPDEL_H_
#define PARSE_POPDEL_H_

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

#include "utils_popdel.h"
#include "popdel_profile/profile_parameter_parsing_popdel.h"
#include "popdel_call/parameter_parsing_popdel_call.h"
#include "popdel_view/popdel_view_parameter_parsing.h"

using namespace seqan;

// -----------------------------------------------------------------------------
// Function printHeader()
// -----------------------------------------------------------------------------
// Print the tool header and description.
void printHeader(const ArgumentParser & parser, int argc, char const ** argv)
{
    std::ostream_iterator<char> out(std::cout);
    CharString name = getName(parser._toolDoc);
    CharString shortDescription = getShortDescription(parser._toolDoc);
    // Print the header.
    std::cout << name;
    if (!empty(shortDescription))
        std::cout << " - " << shortDescription;
    std::cout << "\n";
    unsigned len = length(name) + (empty(shortDescription) ? 0 : 3) + length(shortDescription);
    std::fill_n(out, len, '=');
    std::cout << "\n\n";
    // Print the command.
    std::ostringstream cmd;
    cmd << "Command:  \n\n";
    for (int i = 1; i < argc; ++i)
        cmd << argv[i] << " ";
    cmd << std::endl;
    printStatus(cmd);
}
// -----------------------------------------------------------------------------
// Function parseCommandLine()
// -----------------------------------------------------------------------------
// Initialize the parser the command line or print help.
template<typename TParams>
ArgumentParser::ParseResult parseCommandLine(TParams & params, int argc, char const ** argv)
{
    // Setup the parser.
    ArgumentParser parser("PopDel");
    setupParser(parser, params);
    // Parse the command line and write error messages to error stream.
    std::ostringstream errorStream;
    ArgumentParser::ParseResult res = parse(parser, argc, argv, std::cout, errorStream);
    // Print full help and return if -H option is set, even if another error occurred.
    if (isSet(parser, "full-help") || isSet(parser, "fullHelp"))
    {
        setHiddenOptions(parser, false, params);
        printHelp(parser);
        return ArgumentParser::PARSE_HELP;
    }
    // If error occurred in parsing, now print the error message and return.
    if (res != ArgumentParser::PARSE_OK)
    {
        std::cerr << errorStream.str();
        return res;
    }
    // Print header and command to std::out.
    printHeader(parser, argc, argv);
    // Get the required arguments' and parameters' values.
    getParameterValues(params, parser);
    getArgumentValues(params, parser);
    return ArgumentParser::PARSE_OK;
}
ArgumentParser::ParseResult parseCommandLine(PopDelViewParameters & params, int argc, char const ** argv)
{
    // Setup the parser.
    ArgumentParser parser("PopDel");
    setupParser(parser, params);
    // Parse the command line and write error messages to error stream.
    std::ostringstream errorStream;
    ArgumentParser::ParseResult res = parse(parser, argc, argv, std::cout, errorStream);
    // Print full help and return if -H option is set, even if another error occurred.
    if (isSet(parser, "full-help") || isSet(parser, "fullHelp"))
    {
        setHiddenOptions(parser, false, params);
        printHelp(parser);
        return ArgumentParser::PARSE_HELP;
    }
    // If error occurred in parsing, now print the error message and return.
    if (res != ArgumentParser::PARSE_OK)
    {
        std::cerr << errorStream.str();
        return res;
    }
    // Get the required arguments' and parameters' values.
    getParameterValues(params, parser);
    getArgumentValues(params, parser);
    return ArgumentParser::PARSE_OK;
}
// ==========================================================================
// Function printHelp()
// ==========================================================================
// Print the help lines.
void printHelp(char const * name)
{
    std::cerr << "PopDel - Population-wide deletion calling" << std::endl;
    std::cerr << "=========================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << name << " COMMAND\033[0m [\033[4mOPTIONS\033[0m]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mCOMMANDS\033[0m" << std::endl;
    std::cerr << "    \033[1mprofile\033[0m   Write a profile of read pair distances for a given BAM or CRAM file." << std::endl;
    std::cerr << "    \033[1mcall\033[0m      Call deletions based on read pair distance profiles of multiple samples."
              << std::endl;
    std::cerr << "    \033[1mview\033[0m      View a profile of read pair distances in human readable format." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << "PopDel" << " version: " << VERSION << std::endl;
    std::cerr << "    Last update " << DATE << std::endl;
    std::cerr << "    Contact: Sebastian Niehus (Sebastian.Niehus[at]ukr.de)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try `" << name << " COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
}
inline int checkParser(const ArgumentParser::ParseResult & res)
{
    if (res == ArgumentParser::PARSE_HELP ||
        res == ArgumentParser::PARSE_VERSION ||
        res ==  ArgumentParser::PARSE_WRITE_CTD ||
        res == ArgumentParser::PARSE_EXPORT_HELP)
        return 0;
    else if (res != ArgumentParser::PARSE_OK)
        return 1;
    else
        return -1;
}

#endif /* PARSE_POPDEL_H_ */
