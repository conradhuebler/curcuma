#pragma once

#include "src/core/energycalculator.h"

#include <Eigen/Dense>
#include <vector>

using Vector = Eigen::VectorXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

class LBFGS {
public:
    explicit LBFGS(int m = 10);

    void initialize(int atoms, const Vector& initial_x);

    Vector step();
    Vector LBFGSStep();
    Vector DIISStep();
    Vector RFOStep();

    Vector getCurrentSolution() const;
    double getCurrentObjectiveValue() const;
    double Energy() const { return m_energy; }
    Vector getCurrentGradient() const;
    bool isConverged() { return false; }
    void setEnergyCalculator(EnergyCalculator* interface)
    {
        m_interface = interface;
    }

    void setConstrains(const Vector& vector) { m_constrains = vector; }

    bool isError() { return m_error; }
    void setOptimMethod(int method)
    {
        m_method = method;
    }

    void setHessian(const Matrix& hessian);
    void setLambda(double lambda) { m_lambda = lambda; }
    void setMasses(const std::vector<double>& masses) { m_masses = masses; }
    void setDIIS(int hist, int start)
    {
        m_diis_hist = hist;
        m_diis_start = start;
    }

private:
    double lineSearchRFO(const Vector& x, const Vector& p, double c1 = 1e-4, double c2 = 0.9);
    double line_search_backtracking(const Vector& x, const Vector& p, double alpha_init, double rho, double c1 = 1e-4);

    Vector diisExtrapolation();

    double getEnergyGradient(const Vector& x);
    double getEnergy(const Vector& x);
    void updateHessian();
    Matrix sr1Update(const Eigen::VectorXd& s, const Eigen::VectorXd& y);
    Eigen::MatrixXd update_inverse_hessian_approximation(Eigen::MatrixXd& B, const Eigen::VectorXd& deltaX, const Eigen::VectorXd& deltaG);

    EnergyCalculator* m_interface;

    bool m_error = false;
    double m_lambda = 0.1;
    int m_atoms = 0;
    int m;
    int stepCount;
    int m_method = 1;
    int m_diis_hist = 10, m_diis_start = 10;
    double m_energy, m_dE, m_step_size = 2, m_last_step_size = 2;
    Matrix m_hessian, m_eigenvectors, m_hess_inv;
    Vector x;
    Vector m_gradient;
    Vector p;
    Vector m_constrains, m_eigenvals;
    std::vector<Vector> m_solutions;
    std::vector<Vector> m_errors;

    std::vector<Vector> s_list;
    std::vector<Vector> y_list;
    std::vector<double> rho_list;
    std::vector<double> m_masses;
};
