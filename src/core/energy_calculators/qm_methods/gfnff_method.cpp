/*
 * < Native GFN-FF Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "gfnff_method.h"
#include "src/tools/general.h"

// Minimal stub implementation for GFN-FF method
GFNFFMethod::GFNFFMethod(const json& config)
    : m_gfnff(std::make_unique<GFNFF>(config))
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_parameters = config;
}

bool GFNFFMethod::setMolecule(const Mol& mol) {
    m_molecule = mol;
    m_initialized = true;
    return m_gfnff->QMInterface::InitialiseMolecule(mol);
}

bool GFNFFMethod::updateGeometry(const Matrix& geometry) {
    m_molecule.m_geometry = geometry;
    return m_gfnff->QMInterface::UpdateMolecule(geometry);
}

double GFNFFMethod::calculateEnergy(bool gradient)
{
    m_last_energy = m_gfnff->Calculation(gradient);
    m_calculation_done = true;
    return m_last_energy;
}

Matrix GFNFFMethod::getGradient() const {
    return m_gfnff->Gradient();
}

Vector GFNFFMethod::getCharges() const {
    return m_gfnff->Charges();
}

Vector GFNFFMethod::getBondOrders() const {
    return m_gfnff->BondOrders();
}

Position GFNFFMethod::getDipole() const {
    return Position{0.0, 0.0, 0.0};
}

void GFNFFMethod::setThreadCount(int threads) {
    m_thread_count = threads;
}

void GFNFFMethod::setParameters(const json& params) {
    m_parameters = params;
    if (m_gfnff) {
        m_gfnff->setParameters(params);
    }
}

json GFNFFMethod::getParameters() const {
    return m_parameters;
}

bool GFNFFMethod::hasError() const {
    return m_has_error;
}

void GFNFFMethod::clearError() {
    m_has_error = false;
    m_error_message.clear();
}

std::string GFNFFMethod::getErrorMessage() const {
    return m_error_message;
}