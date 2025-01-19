#include "interface/abstract_interface.h"

#include "GFN.hpp"

#include <string>
#include <vector>

#include "ulyssesinterface.h"

bool UlyssesInterface::InitialiseMolecule(const Molecule& molecule)
{
    reset();

    m_mol = new Molecule(molecule); // Konvertierung zu Ulysses Molecule
    m_initialised = true;

    return setupCalculation();
}

bool UlyssesInterface::InitialiseMolecule(const int* attyp, const double* coord,
    const int natoms, const double charge, const int spin)
{
}

bool UlyssesInterface::setupCalculation()
{
}

double UlyssesInterface::Calculate(double* gradient, bool verbose)
{
    if (!m_calculator_setup)
        return 0.0;

    double energy = m_electron->Calculate(0); // 0 für Standard-Berechnung

    // Lade Charges und Polarizabilities
    m_charges = m_electron->getQAtoms();

    m_electron->AtomicPolarizabilities(m_polarizabilities, m_charges);
    m_electron->TotalPolarizability(m_total_polarizability, m_charges);

    if (gradient != nullptr) {
        // Gradient-Berechnung implementieren wenn nötig
    }

    return energy;
}

/*
std::vector<double> UlyssesInterface::getCharges() const
{
    return m_charges;
}

std::vector<double> UlyssesInterface::getPolarizabilities() const
{
    return m_polarizabilities;
}

double UlyssesInterface::getTotalPolarizability() const
{
    return m_total_polarizability;
}

std::vector<size_t> UlyssesInterface::getAOBasisInfo() const
{
    if(!m_calculator_setup)
        return std::vector<size_t>();

    return m_basis->AtomNAOs(m_mol->Atoms());
}

void UlyssesInterface::reset()
{
    delete m_mol;
    delete m_basis;
    delete m_electron;

    m_mol = nullptr;
    m_basis = nullptr;
    m_electron = nullptr;

    m_initialised = false;
    m_calculator_setup = false;

    m_charges.clear();
    m_polarizabilities.clear();
    m_total_polarizability = 0.0;
}
*/
bool UlyssesInterface::UpdateMolecule(const Molecule& molecule)
{
    // Implementiere Molekül-Update
    return true;
}

bool UlyssesInterface::UpdateMolecule(const double* coord)
{
    // Implementiere Koordinaten-Update
    return true;
}
