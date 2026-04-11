#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#include "src/tools/cli_parser.h"

using json = nlohmann::json;

void testCLI2Json()
{
    // Test case 1: No arguments
    {
        char* argv[] = { (char*)"curcuma" };
        json result = CLIUtils::CLI2Json(1, argv);
        assert(result.empty());
    }

    // Test case 2: Single flag (modern dash handling)
    {
        char* argv[] = { (char*)"curcuma", (char*)"--md", (char*)"--flag" };
        json result = CLIUtils::CLI2Json(3, argv);
        // "md" is mapped to "simplemd"
        assert(result["simplemd"]["flag"] == true);
        assert(result["md"]["flag"] == true);
    }

    // Test case 3: Hyphen to underscore (canonicalization)
    {
        char* argv[] = { (char*)"curcuma", (char*)"-opt", (char*)"--max-iter", (char*)"100" };
        json result = CLIUtils::CLI2Json(4, argv);
        assert(result["opt"]["max_iter"] == 100);
    }

    // Test case 4: Triple dashes and mixed formats
    {
        char* argv[] = { (char*)"curcuma", (char*)"---sp", (char*)"-method", (char*)"uff", (char*)"--threads", (char*)"4" };
        json result = CLIUtils::CLI2Json(6, argv);
        assert(result["opt"]["method"] == "uff");
        assert(result["global"]["method"] == "uff");
        assert(result["global"]["threads"] == 4);
    }

    // Test case 5: Dotted parameters with prefix stripping
    {
        char* argv[] = { (char*)"curcuma", (char*)"-md", (char*)"-md.max_time", (char*)"10.5" };
        json result = CLIUtils::CLI2Json(4, argv);
        // Should be at top level of simplemd, not nested
        assert(result["simplemd"]["max_time"] == 10.5);
    }

    // Test case 6: Dotted parameters routing to other modules
    {
        char* argv[] = { (char*)"curcuma", (char*)"-opt", (char*)"-rmsd.method", (char*)"subspace" };
        json result = CLIUtils::CLI2Json(4, argv);
        assert(result["rmsd"]["method"] == "subspace");
    }

    // Test case 7: Positional arguments handling
    {
        char* argv[] = { (char*)"curcuma", (char*)"-sp", (char*)"molecule.xyz", (char*)"extra_arg" };
        json result = CLIUtils::CLI2Json(4, argv);
        // "sp" is mapped to "opt"
        assert(result["opt"]["input_file"] == "molecule.xyz");
        assert(result["positional"].is_array());
        assert(result["positional"].size() == 2);
        assert(result["positional"][0] == "molecule.xyz");
        assert(result["positional"][1] == "extra_arg");
    }

    std::cout << "All CLI parsing tests passed!" << std::endl;
}

int main()
{
    testCLI2Json();
    return 0;
}
