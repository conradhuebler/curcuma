#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// Data structure to hold a parsed parameter definition
struct Parameter {
    std::string module;
    std::string name;
    std::string type;
    std::string defaultValue;
    std::string helpText;
    std::string category;
    std::string aliases;
};

// Function to trim leading/trailing whitespace and quotes
std::string trim(const std::string& str)
{
    const std::string whitespace = " \t\"";
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return "";
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

// Claude Generated (Sep 2026): tokenizer for PARAM(...) macros, see main() step 2.
struct ParsedMacro {
    std::vector<std::string> args;   ///< top-level arguments, raw text
};

/// Find every PARAM( ... ) in a block of C++ text, skipping strings, char literals and comments.
static void extract_param_macros(const std::string& t, std::vector<ParsedMacro>& out, std::vector<size_t>& starts)
{
    enum State { Code, Str, Chr, LineComment, BlockComment };
    State st = Code;
    const size_t n = t.size();
    for (size_t i = 0; i < n; ++i) {
        const char c = t[i];
        if (st == LineComment) { if (c == '\n') st = Code; continue; }
        if (st == BlockComment) { if (c == '*' && i + 1 < n && t[i + 1] == '/') { st = Code; ++i; } continue; }
        if (st == Str) { if (c == '\\') ++i; else if (c == '"') st = Code; continue; }
        if (st == Chr) { if (c == '\\') ++i; else if (c == '\'') st = Code; continue; }
        if (c == '/' && i + 1 < n && t[i + 1] == '/') { st = LineComment; ++i; continue; }
        if (c == '/' && i + 1 < n && t[i + 1] == '*') { st = BlockComment; ++i; continue; }
        if (c == '"') { st = Str; continue; }
        if (c == '\'') { st = Chr; continue; }
        // identifier PARAM at a word boundary, followed by optional whitespace and '('
        if (c != 'P' || t.compare(i, 5, "PARAM") != 0) continue;
        if (i > 0 && (std::isalnum(static_cast<unsigned char>(t[i - 1])) || t[i - 1] == '_')) continue;
        size_t k = i + 5;
        if (k < n && (std::isalnum(static_cast<unsigned char>(t[k])) || t[k] == '_')) continue;
        while (k < n && std::isspace(static_cast<unsigned char>(t[k]))) ++k;
        if (k >= n || t[k] != '(') continue;
        // Collect the arguments up to the matching ')'.
        ParsedMacro m;
        std::string cur;
        int depth = 0;           // nesting of (), {}, [] inside the macro
        State s2 = Code;
        size_t e = k + 1;
        bool closed = false;
        for (; e < n; ++e) {
            const char d = t[e];
            if (s2 == LineComment) { if (d == '\n') s2 = Code; continue; }
            if (s2 == BlockComment) { if (d == '*' && e + 1 < n && t[e + 1] == '/') { s2 = Code; ++e; } continue; }
            if (s2 == Str || s2 == Chr) {
                cur += d;
                if (d == '\\' && e + 1 < n) { cur += t[++e]; continue; }
                if ((s2 == Str && d == '"') || (s2 == Chr && d == '\'')) s2 = Code;
                continue;
            }
            if (d == '/' && e + 1 < n && t[e + 1] == '/') { s2 = LineComment; ++e; continue; }
            if (d == '/' && e + 1 < n && t[e + 1] == '*') { s2 = BlockComment; ++e; continue; }
            if (d == '"') { s2 = Str; cur += d; continue; }
            if (d == '\'') { s2 = Chr; cur += d; continue; }
            if (d == '(' || d == '{' || d == '[') ++depth;
            if (d == ')' && depth == 0) { m.args.push_back(cur); closed = true; break; }
            if (d == ')' || d == '}' || d == ']') --depth;
            if (d == ',' && depth == 0) { m.args.push_back(cur); cur.clear(); continue; }
            cur += d;
        }
        starts.push_back(i);
        if (!closed) m.args.clear();   // unterminated: reported as malformed
        out.push_back(m);
        i = closed ? e : n;
    }
}

/// Join the string literal(s) of one argument ("a" "b" -> a b, escapes kept verbatim).
/// Returns false when the argument contains anything but string literals and whitespace.
static bool join_string_literals(const std::string& a, std::string& out)
{
    out.clear();
    bool any = false;
    for (size_t i = 0; i < a.size(); ++i) {
        const char c = a[i];
        if (std::isspace(static_cast<unsigned char>(c))) continue;
        if (c != '"') return false;
        any = true;
        for (++i; i < a.size() && a[i] != '"'; ++i) {
            if (a[i] == '\\' && i + 1 < a.size()) out += a[i++];
            out += a[i];
        }
        if (i >= a.size()) return false;
    }
    return any;
}

/// Turn the six macro arguments into a Parameter, with the same field conventions as the old
/// regex parser (trimmed name/type/default/category; help = literal content; aliases = raw text
/// between the braces, untrimmed).
static bool macro_to_parameter(const ParsedMacro& m, const std::string& module, Parameter& p)
{
    if (m.args.size() != 6) return false;
    std::string help, cat;
    if (!join_string_literals(m.args[3], help) || !join_string_literals(m.args[4], cat)) return false;
    const std::string& al = m.args[5];
    const size_t ob = al.find('{'), cb = al.rfind('}');
    if (ob == std::string::npos || cb == std::string::npos || cb < ob) return false;
    p.module = module;
    p.name = trim(m.args[0]);
    p.type = trim(m.args[1]);
    p.defaultValue = trim(m.args[2]);
    p.helpText = trim(help);
    p.category = trim(cat);
    p.aliases = al.substr(ob + 1, cb - ob - 1);
    return !p.name.empty() && !p.type.empty();
}

// Function to generate the C++ header file content - Claude Generated (fixed)
std::string generate_header_content(const std::vector<Parameter>& params)
{
    std::stringstream ss;
    ss << "// THIS FILE IS AUTO-GENERATED BY curcuma_param_parser. DO NOT EDIT MANUALLY.\n\n";
    ss << "#pragma once\n";
    ss << "#include \"src/core/parameter_registry.h\"\n";
    ss << "#include <string>\n\n";
    ss << "// This function will be called to populate the registry\n";
    ss << "inline void initialize_generated_registry() {\n";
    ss << "    auto& registry = ParameterRegistry::getInstance();\n\n";

    for (const auto& p : params) {
        ss << "    registry.addDefinition(\"" << p.module << "\", {\n";
        ss << "        \"" << p.name << "\",\n";
        ss << "        \"" << p.module << "\",\n";
        ss << "        ParamType::" << p.type << ",\n";

        // Type-appropriate default value formatting - Claude Generated
        if (p.type == "String") {
            ss << "        std::string(\"" << p.defaultValue << "\"),\n";
        } else if (p.type == "Bool") {
            ss << "        " << p.defaultValue << ",\n";
        } else {
            // Int or Double - no quotes needed
            ss << "        " << p.defaultValue << ",\n";
        }

        ss << "        \"" << p.helpText << "\",\n";
        ss << "        \"" << p.category << "\",\n";

        // Handle aliases - Claude Generated (Fixed October 2025 - support multiple aliases)
        if (p.aliases.empty()) {
            ss << "        {}\n";
        } else {
            // DON'T trim quotes here - we need them for parsing!
            // Only trim whitespace
            std::string aliases_str = p.aliases;
            // Remove leading/trailing whitespace only (not quotes)
            size_t start = aliases_str.find_first_not_of(" \t");
            size_t end = aliases_str.find_last_not_of(" \t");
            if (start == std::string::npos) {
                ss << "        {}\n";
            } else {
                aliases_str = aliases_str.substr(start, end - start + 1);

                // Parse multiple comma-separated alias strings
                std::vector<std::string> alias_list;
                std::string current_alias;
                bool in_quotes = false;

                for (char c : aliases_str) {
                    if (c == '"') {
                        in_quotes = !in_quotes;
                        // Don't add quotes to the alias string itself
                    } else if (c == ',' && !in_quotes) {
                        // End of one alias, start of next
                        std::string trimmed_alias = trim(current_alias);
                        if (!trimmed_alias.empty()) {
                            alias_list.push_back(trimmed_alias);
                        }
                        current_alias.clear();
                    } else if (in_quotes || (c != ' ' && c != '\t')) {
                        // Add character if inside quotes OR if it's not whitespace outside quotes
                        current_alias += c;
                    }
                }

                // Add last alias
                std::string trimmed_alias = trim(current_alias);
                if (!trimmed_alias.empty()) {
                    alias_list.push_back(trimmed_alias);
                }

                // Generate initializer list
                if (alias_list.empty()) {
                    ss << "        {}\n";
                } else {
                    ss << "        {";
                    for (size_t i = 0; i < alias_list.size(); ++i) {
                        ss << "\"" << alias_list[i] << "\"";
                        if (i + 1 < alias_list.size()) {
                            ss << ", ";
                        }
                    }
                    ss << "}\n";
                }
            }
        }

        ss << "    });\n";
    }

    ss << "}\n";
    return ss.str();
}

int main(int argc, char* argv[])
{
    std::vector<std::string> input_files;
    std::string output_file;

    // 1. Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--inputs" && i + 1 < argc) {
            // Collect all subsequent arguments as input files until another flag is found
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                input_files.push_back(argv[++i]);
            }
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        }
    }

    if (output_file.empty()) {
        std::cerr << "Error: --output file path not specified." << std::endl;
        return 1;
    }
    if (input_files.empty()) {
        std::cerr << "Warning: No --inputs specified." << std::endl;
    }

    // 2. Parse all input files
    //
    // Claude Generated (Sep 2026, docs/MULTI_GPU_GAPS.md X-1): the PARAM macros are found by a
    // small tokenizer instead of a line-accumulating regex. The regex failed on multi-line PARAMs
    // whose help text is several adjacent string literals ("a" "b") or contains ')' - and a
    // failed match left the rest of that PARAM in the buffer, so the following PARAMs of the block
    // were dropped too (26 names in eeq_solver.h and gfnff.h, silently absent from -help,
    // flat-flag routing and -export_run). The tokenizer tracks strings, character literals and
    // comments, takes the macro up to its matching ')', splits the six top-level arguments and
    // joins adjacent string literals. A PARAM inside a comment is not a definition.
    std::vector<Parameter> all_params;
    std::regex begin_regex(R"(BEGIN_PARAMETER_DEFINITION\s*\(\s*(\w+)\s*\))");
    std::regex end_regex(R"(END_PARAMETER_DEFINITION)");
    int malformed = 0;

    for (const auto& filepath : input_files) {
        // Claude Generated (October 2025): Skip parameter_macros.h - it contains macro definitions only, not parameter declarations
        if (filepath.find("parameter_macros.h") != std::string::npos) {
            continue;
        }

        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Warning: Could not open file " << filepath << "\n";
            continue;
        }

        // Collect the text of every BEGIN_PARAMETER_DEFINITION ... END_PARAMETER_DEFINITION block,
        // with the line number of each block line for the diagnostics.
        std::string line;
        int line_number = 0;
        bool in_block = false;
        std::string module, block;
        std::vector<int> block_line_of_char;
        auto flush_block = [&]() {
            std::vector<ParsedMacro> macros;
            std::vector<size_t> starts;
            extract_param_macros(block, macros, starts);
            for (size_t k = 0; k < macros.size(); ++k) {
                const int at = starts[k] < block_line_of_char.size() ? block_line_of_char[starts[k]] : line_number;
                Parameter prm;
                if (!macro_to_parameter(macros[k], module, prm)) {
                    std::cerr << "Warning: Malformed PARAM in " << filepath << " around line " << at << "\n";
                    ++malformed;
                    continue;
                }
                all_params.push_back(prm);
            }
            block.clear();
            block_line_of_char.clear();
        };
        while (std::getline(file, line)) {
            ++line_number;
            std::smatch match;
            if (!in_block && std::regex_search(line, match, begin_regex)) {
                in_block = true;
                module = match[1];
                continue;
            }
            if (in_block && std::regex_search(line, match, end_regex)) {
                flush_block();
                in_block = false;
                module.clear();
                continue;
            }
            if (in_block) {
                block += line;
                block += '\n';
                block_line_of_char.resize(block.size(), line_number);
            }
        }
        if (in_block)
            flush_block();
    }
    if (malformed > 0)
        std::cerr << "Warning: " << malformed << " malformed PARAM definition(s) were skipped\n";

    // 3. Generate the header file content
    std::string header_content = generate_header_content(all_params);

    // 4. Write the output file
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open output file " << output_file << " for writing." << std::endl;
        return 1;
    }
    out << header_content;
    out.close();

    std::cout << "Successfully generated " << output_file << " with " << all_params.size() << " parameter definitions." << std::endl;

    return 0;
}
