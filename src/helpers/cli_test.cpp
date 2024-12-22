#include <cassert>
#include <iostream>

#include "src/tools/general.h"

#include "json.hpp" // Include the JSON library
#include "src/core/global.h" // Include the header file with the CLI2Json function

using json = nlohmann::json;

void testCLI2Json()
{
    // Test case 1: No arguments
    {
        char* argv[] = { (char*)"program" };
        json result = CLI2Json(1, argv);
        assert(result.empty());
    }

    // Test case 2: Single flag
    {
        char* argv[] = { (char*)"program", (char*)"-keyword", (char*)"-flag" };
        json result = CLI2Json(3, argv);
        assert(result["keyword"]["flag"] == true);
    }

    // Test case 3: Flag with value
    {
        char* argv[] = { (char*)"program", (char*)"-keyword", (char*)"-flag", (char*)"value" };
        json result = CLI2Json(4, argv);
        assert(result["keyword"]["flag"] == "value");
    }

    // Test case 4: Flag with true
    {
        char* argv[] = { (char*)"program", (char*)"-keyword", (char*)"-flag", (char*)"true" };
        json result = CLI2Json(4, argv);
        assert(result["keyword"]["flag"] == true);
    }

    // Test case 5: Flag with false
    {
        char* argv[] = { (char*)"program", (char*)"-keyword", (char*)"-flag", (char*)"false" };
        json result = CLI2Json(4, argv);
        assert(result["keyword"]["flag"] == false);
    }

    // Test case 6: Flag with number
    {
        char* argv[] = { (char*)"program", (char*)"-keyword", (char*)"-flag", (char*)"123.45" };
        json result = CLI2Json(4, argv);
        assert(result["keyword"]["flag"] == 123.45);
    }

    // Test case 7: Flag with vector
    {
        char* argv[] = { (char*)"program", (char*)"-keyword", (char*)"-flag", (char*)"1,2,3" };
        json result = CLI2Json(4, argv);
        assert(result["keyword"]["flag"] == "1,2,3");
    }

    std::cout << "All tests passed!" << std::endl;
}

int main()
{
    testCLI2Json();
    return 0;
}