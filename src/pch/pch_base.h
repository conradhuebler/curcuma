/*
 * Base Precompiled Header - Standard Library + External Dependencies
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "src/global_config.h"

// Standard C++ Library Headers (most frequently used)
#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

// C++17 specific headers
#ifdef C17
#include <filesystem>
#endif

// Third-party libraries that are safe to precompile
#include "json.hpp"
#include <Eigen/Dense>
#include <fmt/core.h>