/*
 * <Unified XYZ Comment Parser Implementation>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Unifies 10 duplicate setXYZComment_X functions
 */

#include "xyz_comment_parser.h"
#include "molecule.h"
#include "src/tools/general.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace curcuma {

// ============================================================================
// Public Interface
// ============================================================================

bool XYZCommentParser::parseComment(const std::string& comment, Molecule& mol)
{
    if (isWhitespaceOnly(comment)) {
        return true; // Empty comment is valid, just ignore
    }

    FormatType format = detectFormat(comment);
    std::vector<std::string> tokens = tokenize(comment);

    // Debug: Show that new parser is being used (remove after testing)
    if (comment.find("water") != std::string::npos || comment.find("Energy") != std::string::npos) {
        std::cout << "[DEBUG] XYZCommentParser: Processing \"" << comment << "\" as format " << static_cast<int>(format) << std::endl;
    }

    switch (format) {
    case FormatType::CURCUMA_NATIVE:
        return parseCurcumaNative(tokens, mol);

    case FormatType::XTB_GNORM:
        return parseXTBGnorm(tokens, mol);

    case FormatType::GAUSSIAN_SCF:
        return parseGaussianSCF(tokens, mol);

    case FormatType::ORCA_ENERGY:
        return parseORCAEnergy(comment, mol);

    case FormatType::SIMPLE_ENERGY:
        return parseSimpleEnergy(tokens, mol);

    case FormatType::GENERIC_6_FIELD:
        return parseGeneric6Field(tokens, mol);

    case FormatType::GENERIC_7_FIELD:
        return parseGeneric7Field(tokens, mol);

    case FormatType::GENERIC_MULTI:
        return parseGenericMulti(tokens, mol);

    case FormatType::EMPTY:
        return true; // Valid but no-op

    case FormatType::UNKNOWN:
    default:
        // Fallback: try to extract any number as energy
        return parseGenericMulti(tokens, mol);
    }
}

// ============================================================================
// Format Detection
// ============================================================================

XYZCommentParser::FormatType XYZCommentParser::detectFormat(const std::string& comment)
{
    if (isWhitespaceOnly(comment)) {
        return FormatType::EMPTY;
    }

    std::vector<std::string> tokens = tokenize(comment);

    // Curcuma native format: "** Energy = X ** Charge = Y ** Curcuma"
    if (comment.find("Curcuma") != std::string::npos && tokens.size() >= 8) {
        return FormatType::CURCUMA_NATIVE;
    }

    // XTB format: " energy: X gnorm: Y xtb: Z" (exactly 7 tokens with "gnorm:")
    if (tokens.size() == 7 && tokens.size() > 2 && tokens[2] == "gnorm:") {
        return FormatType::XTB_GNORM;
    }

    // Gaussian SCF format: "SCF done X ..." (exactly 4 tokens starting with "SCF done")
    if (tokens.size() == 4 && tokens.size() >= 2 && tokens[0] == "SCF" && tokens[1] == "done") {
        return FormatType::GAUSSIAN_SCF;
    }

    // ORCA format: Contains specific signature
    if (comment.find("Coordinates from ORCA-job input E") != std::string::npos) {
        return FormatType::ORCA_ENERGY;
    }

    // Simple energy: Single token that is a number
    if (tokens.size() == 1 && isNumber(tokens[0])) {
        return FormatType::SIMPLE_ENERGY;
    }

    // Generic formats based on token count (preserve legacy dispatcher behavior)
    if (tokens.size() == 6) {
        return FormatType::GENERIC_6_FIELD;
    }

    if (tokens.size() == 7) {
        // 7 fields but NOT XTB format (gnorm already handled above)
        return FormatType::GENERIC_7_FIELD;
    }

    // Multi-field fallback (extract first number found)
    if (tokens.size() > 1) {
        return FormatType::GENERIC_MULTI;
    }

    return FormatType::UNKNOWN;
}

// ============================================================================
// Specific Format Parsers (preserve exact legacy behavior)
// ============================================================================

bool XYZCommentParser::parseCurcumaNative(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Legacy behavior: Try two different energy positions
    try {
        mol.setEnergy(std::stod(tokens[4]));
        mol.setCharge(std::stod(tokens[9]));
    } catch (const std::invalid_argument& what) {
        try {
            mol.setEnergy(std::stod(tokens[3]));
            mol.setCharge(std::stod(tokens[8]));
        } catch (const std::invalid_argument& what) {
            // Ignore parsing errors, keep existing values
        }
    }
    return true;
}

bool XYZCommentParser::parseXTBGnorm(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Legacy setXYZComment_7 with gnorm detection
    // Format: " energy: X gnorm: Y xtb: Z" → extract energy from position [1]
    try {
        mol.setEnergy(extractNumber(tokens[1]));
    } catch (const std::exception&) {
        mol.setEnergy(0);
    }
    return true;
}

bool XYZCommentParser::parseGaussianSCF(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Legacy setXYZComment_4 SCF parsing
    // Format: "SCF done <energy> <other>" → extract energy from position [2]
    if (tokens[0] == "SCF" && tokens[1] == "done") {
        try {
            mol.setEnergy(extractNumber(tokens[2]));
        } catch (const std::exception&) {
            mol.setEnergy(0);
        }
    } else {
        // Alternative 4-field format: name at [0], energy at [3]
        mol.setName(tokens[0]);
        if (tokens[3].empty()) {
            mol.setEnergy(0);
        } else {
            try {
                mol.setEnergy(extractNumber(tokens[3]));
            } catch (const std::exception&) {
                mol.setEnergy(0);
            }
        }
    }
    return true;
}

bool XYZCommentParser::parseORCAEnergy(const std::string& comment, Molecule& mol)
{
    // Extract energy after "E " in ORCA comment
    size_t pos = comment.find(" E ");
    if (pos != std::string::npos) {
        std::string energy_part = comment.substr(pos + 3);
        std::vector<std::string> tokens = tokenize(energy_part);
        if (!tokens.empty()) {
            try {
                mol.setEnergy(extractNumber(tokens[0]));
            } catch (const std::exception&) {
                mol.setEnergy(0);
            }
        }
    }
    return true;
}

bool XYZCommentParser::parseSimpleEnergy(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Single numeric field - set as energy directly
    try {
        mol.setEnergy(extractNumber(tokens[0]));
    } catch (const std::exception&) {
        mol.setEnergy(0);
    }
    return true;
}

bool XYZCommentParser::parseGeneric6Field(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Legacy setXYZComment_6: energy at position [3]
    try {
        mol.setEnergy(extractNumber(tokens[3]));
    } catch (const std::exception&) {
        mol.setEnergy(0);
    }
    return true;
}

bool XYZCommentParser::parseGeneric7Field(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Legacy setXYZComment_7 non-gnorm path: energy at position [4]
    try {
        mol.setEnergy(extractNumber(tokens[4]));
    } catch (const std::exception&) {
        mol.setEnergy(0);
    }
    return true;
}

bool XYZCommentParser::parseGenericMulti(const std::vector<std::string>& tokens, Molecule& mol)
{
    // Legacy setXYZComment_0/1/2/3 behavior: find first valid number
    for (const std::string& token : tokens) {
        if (isNumber(token)) {
            try {
                mol.setEnergy(extractNumber(token));
                break; // Take first valid number found
            } catch (const std::exception&) {
                continue; // Try next token
            }
        }
    }
    return true;
}

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<std::string> XYZCommentParser::tokenize(const std::string& comment)
{
    // Use existing SplitString function for consistency
    return Tools::SplitString(comment);
}

bool XYZCommentParser::isNumber(const std::string& str)
{
    // Use existing Tools::isDouble for consistency
    return Tools::isDouble(str);
}

bool XYZCommentParser::isWhitespaceOnly(const std::string& str)
{
    return str.empty() || std::all_of(str.begin(), str.end(), [](unsigned char c) { return std::isspace(c); });
}

double XYZCommentParser::extractNumber(const std::string& str, double fallback)
{
    try {
        return std::stod(str);
    } catch (const std::exception&) {
        return fallback;
    }
}

} // namespace curcuma