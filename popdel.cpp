#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

#include "parse_popdel.h"
#include "workflow_popdel.h"
#include "popdel_profile/bam_window_iterator_popdel.h"
#include "popdel_profile/parameter_estimation_popdel.h"
#include "popdel_profile/profile_parameter_parsing_popdel.h"

using namespace seqan;

// ==========================================================================
// Function main()
// ==========================================================================
int main(int argc, char const ** argv)
{
    const char * prog_name = argv[0];
    if (argc < 2)
    {
        printHelp(prog_name);
        return 1;
    }
    // Concatenate program name and command.
    const char * command = argv[1];
    std::ostringstream name;
    name << argv[0] << " " << command;
    argv[1] = toCString(name.str());
    ++argv;
    --argc;
    // Execute the specified command.
    if (strcmp(command, "profile") == 0)
    {
        return popdel_profile(argc, argv);
    }
    else if (strcmp(command, "call") == 0)
    {
        return popdel_call(argc, argv);
    }
    else if (strcmp(command, "view") == 0)
    {
        return popdel_view(argc, argv);
    }
    else if (strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0)
    {
        printHelp(prog_name);
        return 0;
    }
    else
    {
        std::cerr << "ERROR: Unknown command: " << command << std::endl;
        printHelp(prog_name);
        return 1;
    }
}
