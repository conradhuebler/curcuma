#pragma once

#include "interface/abstract_interface.h"

#include "GFN.hpp"

#include <string>
#include <vector>

class UlyssesInterface : public QMInterface {
public:
    UlyssesInterface() = default;
    ~UlyssesInterface() = default;

    bool InitialiseMolecule(const Mol& mol);
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin);
    bool UpdateMolecule(const Molecule& molecule);
    bool UpdateMolecule(const double* coord);

    double Calculate(double* gradient = 0, bool verbose = false);

    // Setter Methoden
    void setElectronTemp(double temp) { m_electron_temp = temp; }
    void setSolvent(const std::string& solvent) { m_solvent = solvent; }

    // Getter Methoden
    /*
    virtual Vector Charges() const { return Vector{}; }
    virtual Vector Dipole() const { return Vector{}; }

    virtual Vector BondOrders() const  { return Vector{}; }
    virtual Vector OrbitalEnergies() const { return Vector{}; }
    virtual Vector OrbitalOccupations() const  { return Vector{}; }

    virtual void setMethod(const std::string& method) { m_method = method; }
    */
private:
    void reset();
    bool setupCalculation();

    Mol m_mol;
    Molecule* m_molecule = nullptr;
    BSet* m_basis = nullptr;
    GFN2* m_electron = nullptr;

    bool m_initialised = false;
    bool m_calculator_setup = false;

    double m_electron_temp = 300.0;
    std::string m_solvent = "none";

    std::vector<double> m_charges;
    std::vector<double> m_polarizabilities;
    double m_total_polarizability = 0.0;
};
