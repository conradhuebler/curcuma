/*
 * <Curcuma Logging System Implementation>
 * Copyright (C) 2025 Claude AI - Generated Code
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "src/core/global.h"

// Static member definitions - Claude Generated
int CurcumaLogger::m_verbosity = 1;
bool CurcumaLogger::m_use_colors = true;
CurcumaLogger::OutputFormat CurcumaLogger::m_format = CurcumaLogger::OutputFormat::TERMINAL;
std::chrono::high_resolution_clock::time_point CurcumaLogger::m_start_time;