#include "lbfgs.h"
#include <iostream>

#include "src/core/energycalculator.h"

LBFGS::LBFGS(int m)
    : m(m)
    , stepCount(0)
    , m_energy(0.0)
{
    m_step_size = 2;
}

void LBFGS::initialize(int natoms, const Vector& initial_x)
{
    m_atoms = natoms;
    m_gradient = Eigen::VectorXd::Zero(3 * m_atoms);
    m_constrains = Eigen::VectorXd::Ones(3 * m_atoms);
    m_hessian = Eigen::MatrixXd::Identity(3 * m_atoms, 3 * m_atoms);
    m_hess_inv = Eigen::MatrixXd::Identity(3 * m_atoms, 3 * m_atoms);
    x = initial_x;
    stepCount = 0;

    s_list.clear();
    y_list.clear();
    rho_list.clear();
}

Vector LBFGS::step()
{
    Vector vector = x;
    m_last_step_size = m_step_size;
    if (m_method == 1 || (m_method == 2 && stepCount < m_diis_start)) // LBFGS
        vector = LBFGSStep();
    else if (m_method == 2) // DIIS-Step
    {
        vector = DIISStep();
    } else if (m_method == 3) // RFO-Step
    {
        vector = RFOStep();
    }
    stepCount++;
    if (m_solutions.size() == m_diis_hist) {
        m_solutions.erase(m_solutions.begin());
        m_errors.erase(m_errors.begin());
    }
    return vector;
}
Eigen::MatrixXd LBFGS::update_inverse_hessian_approximation(Eigen::MatrixXd& B, const Eigen::VectorXd& deltaX, const Eigen::VectorXd& deltaG)
{
    double rho = 1.0 / deltaG.dot(deltaX);
    Eigen::VectorXd B_deltaX = B * deltaX;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(B.rows(), B.cols());

    B = (I - rho * deltaX * deltaG.transpose()) * B * (I - rho * deltaG * deltaX.transpose()) + rho * (deltaX * deltaX.transpose());

    return B;
}

void LBFGS::setHessian(const Matrix& hessian)
{
    m_hessian = hessian;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(m_hessian);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    }

    m_eigenvals = solver.eigenvalues();
    m_eigenvectors = solver.eigenvectors();
    m_hess_inv = m_hessian.inverse();
}

double LBFGS::lineSearchRFO(const Vector& x, const Vector& p, double c1, double c2)
{
    double alpha = m_step_size;
    double alpha_prev = m_last_step_size;

    double f_x = m_energy;
    double gdotp = m_gradient.dot(p);
    int maxiter = 100;
    int iter = 0;
    while (std::abs(alpha - alpha_prev) > 1e-2) {
        Vector x_new = x + alpha * p;

        double f_x_new = getEnergy(x_new);

        if (f_x_new > f_x + c1 * alpha * gdotp || (alpha > 1e-10 && f_x_new >= getEnergy(x + alpha_prev * p))) {
            double alpha_temp = alpha;
            alpha = (alpha_prev + alpha) / 2.0;
            alpha_prev = alpha_temp;
        } else {
            getEnergyGradient(x_new);
            Vector grad_new = m_gradient;

            if (std::abs(grad_new.dot(p)) <= -c2 * gdotp) {
                return alpha;
            } else if (grad_new.dot(p) >= 0) {
                double alpha_temp = alpha;
                alpha = (alpha_prev + alpha) / 2.0;
                alpha_prev = alpha_temp;
            } else {
                alpha_prev = alpha;
                alpha *= 2.0;
            }
        }
        iter++;
    }
    return alpha;
}

double LBFGS::line_search_backtracking(const Vector& x, const Vector& p, double alpha_init, double rho, double c1)
{
    double alpha = alpha_init;
    double f_x = getEnergy(x);
    double gradient_dot_p = m_gradient.dot(p);

    while (getEnergy(x + alpha * p) > f_x + c1 * alpha * gradient_dot_p) {
        alpha *= rho;
    }
    return alpha;
}

Vector LBFGS::diisExtrapolation()
{
    int numErrors = m_errors.size();
    Matrix B = Matrix::Zero(numErrors + 1, numErrors + 1);
    Matrix A = Matrix::Zero(numErrors, numErrors);
    Vector rhs = Vector::Zero(numErrors + 1);
    rhs(numErrors) = 1.0;

    for (int i = 0; i < numErrors; ++i) {
        for (int j = 0; j < numErrors; ++j) {
            B(i, j) = (m_hess_inv * m_errors[i]).dot(m_hess_inv * m_errors[j]);
            A(i, j) = B(i, j);
        }
        B(numErrors, i) = 1.0;
        B(i, numErrors) = 1.0;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    double conditionNumber = svd.singularValues()(0) / svd.singularValues().tail(1)(0);

    if (conditionNumber > 1e10) {
        m_solutions.erase(m_solutions.begin());
        m_errors.erase(m_errors.begin());
        std::cout << m_errors.size() << " " << conditionNumber << std::endl;
        int numErrors = m_errors.size();
        Matrix B = Matrix::Zero(numErrors + 1, numErrors + 1);
        Vector rhs = Vector::Zero(numErrors + 1);
        rhs(numErrors) = -1.0;

        for (int i = 0; i < numErrors; ++i) {
            for (int j = 0; j < numErrors; ++j) {
                B(i, j) = (m_hess_inv * m_errors[i]).dot(m_hess_inv * m_errors[j]);
            }
            B(numErrors, i) = -1.0;
            B(i, numErrors) = -1.0;
        }
        Vector coeff = B.colPivHouseholderQr().solve(rhs);

        Vector diisSolution = Vector::Zero(m_solutions[0].size());
        for (int i = 0; i < numErrors; ++i) {
            diisSolution += coeff[i] * m_solutions[i];
        }

        return diisSolution;
    } else {

        Vector coeff = B.colPivHouseholderQr().solve(rhs);

        Vector diisSolution = Vector::Zero(m_solutions[0].size());
        for (int i = 0; i < numErrors; ++i) {
            diisSolution += coeff[i] * m_solutions[i];
        }

        return diisSolution;
    }
}
Vector LBFGS::DIISStep()
{
    double energy = m_energy;

    if (m_solutions.size() == m_diis_hist) {
        m_solutions.erase(m_solutions.begin());
        m_errors.erase(m_errors.begin());
    }
    getEnergyGradient(x);

    m_solutions.push_back(x);
    m_errors.push_back(m_gradient);

    auto tmp = diisExtrapolation();
    x = tmp;
    return x;
}

Vector LBFGS::LBFGSStep()
{
    double energy = m_energy;
    if (stepCount == 0) {
        getEnergyGradient(x);
        p = -m_gradient;
    } else {
        Vector q = m_gradient;
        std::vector<double> alpha_list;

        for (int i = static_cast<int>(s_list.size()) - 1; i >= 0; --i) {
            double alpha = rho_list[i] * s_list[i].dot(q);
            alpha_list.push_back(alpha);
            q -= alpha * y_list[i];
        }

        Vector z = q;

        for (size_t i = 0; i < s_list.size(); ++i) {
            double beta = rho_list[i] * y_list[i].dot(z);
            z += s_list[i] * (alpha_list[s_list.size() - i - 1] - beta);
        }

        p = -z;
    }
    // m_step_size = 2;
    // m_step_size = line_search_backtracking(x,p, m_step_size, 1e-1, 1e-1);
    // m_step_size = lineSearchRFO(x, p);

    Vector x_old = x;
    Vector gradient_old = m_gradient;
    x += m_step_size * p;
    getEnergyGradient(x);
    if ((energy - m_energy) < 0) {
        m_step_size *= 0.9;
        //    std::cout << m_step_size << std::endl;
    } else {
        //  m_step_size *= 1.01;
    }
    Vector s = x - x_old;
    Vector y = m_gradient - gradient_old;

    double rho = 1.0 / y.dot(s);

    m_solutions.push_back(x);
    m_errors.push_back(m_gradient);

    if (s_list.size() == m) {
        s_list.erase(s_list.begin());
        y_list.erase(y_list.begin());
        rho_list.erase(rho_list.begin());
    }

    s_list.push_back(s);
    y_list.push_back(y);
    rho_list.push_back(rho);

    return x;
}

void LBFGS::updateHessian()
{
    Eigen::SelfAdjointEigenSolver<Matrix> solver(m_hessian);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    }

    // Erhalte die Eigenwerte und Eigenvektoren
    m_eigenvals = solver.eigenvalues();
    m_eigenvectors = solver.eigenvectors();

    // Bereite die RFO-Koeffizienten vor
    Vector tau(m_atoms);
    for (int i = 0; i < m_atoms; ++i) {
        tau(i) = 1.0 / (m_eigenvals(i) + m_lambda); // Rationale Funktionsannäherung
    }
    // Aktualisiere das Hessian mit RFO
    Matrix B = Eigen::MatrixXd::Identity(3 * m_atoms, 3 * m_atoms);
    for (int i = 0; i < 3 * m_atoms; i++) {
        B += tau(i) * m_eigenvectors.col(i) * m_eigenvectors.col(i).transpose();
    }
    m_hessian = (m_hessian + B) / 2;
}

Vector LBFGS::RFOStep()
{
    // Wähle den minimalsten/niedrigsten (negativsten) Eigenwert & Eigenvektor
    int minIndex = 0;
    m_eigenvals.minCoeff(&minIndex);
    //  std::cout << minIndex << "  " << m_eigenvals.minCoeff(&minIndex) << std::endl;
    Vector direction = m_eigenvectors.col(minIndex);
    Vector gradient_mass_weighted = m_gradient;
    for (int i = 0; i < m_gradient.size(); ++i) {
        gradient_mass_weighted[i] /= std::sqrt(m_masses[i]);
        //  std::cout << std::sqrt(m_masses[i]) << " "  << std::endl;// 3 Dimensionen (x, y, z) pro Atom
    }

    Vector gradient_normal = m_eigenvectors.transpose() * gradient_mass_weighted;
    // std::cout << gradient_normal.transpose() << std::endl << m_gradient.transpose() << std::endl;
    //  Berechne den Schritt: nutze direction für EVF
    double lambda = m_eigenvals[minIndex]; // kleinster Eigenwert
    double tau = 0.5; // Dämpfungsfaktor zur vorgeschlagenen Schrittgröße

    Vector p_rfo = -tau * (gradient_normal - lambda * direction); // RFO Schritte mit EVF

    // Maßnahme: Line-Search (falls notwendig)
    // if(std::abs(m_step_size  + 1)<1e-10)
    // m_step_size = line_search_backtracking(x,p_rfo, 1, 1e-4);
    //   m_step_size = 1e-4;
    //   std::cout << m_step_size << " ";
    m_step_size = -0.1; // lineSearchRFO(x, p_rfo);
    // m_step_size = 1;
    //    std::cout << m_step_size << " " <<std::endl;

    // Aktualisiere die aktuelle Position
    Vector x_new = x + m_step_size * p_rfo;
    Vector gradient_prev = m_gradient;
    getEnergyGradient(x_new);
    // updateHessian();

    // Update des Hessians (SR1)
    Eigen::VectorXd s = x_new - x;
    Eigen::VectorXd y = m_gradient - gradient_prev;
    Matrix H = sr1Update(s, y);
    setHessian(H);
    // std::cout << ( x - x_new).transpose() << std::endl;
    x = x_new;
    return x_new;
}
Matrix LBFGS::sr1Update(const Eigen::VectorXd& s, const Eigen::VectorXd& y)
{
    Eigen::VectorXd diff = y - m_hessian * s;
    double denom = diff.dot(s);
    if (std::abs(denom) > 1e-10) {
        return m_hessian + m_lambda * (diff * diff.transpose()) / denom;
    }
    return m_hessian;
}
Vector LBFGS::getCurrentSolution() const
{
    return x;
}

double LBFGS::getCurrentObjectiveValue() const
{
    return m_energy;
}

Vector LBFGS::getCurrentGradient() const
{
    return m_gradient;
}

double LBFGS::getEnergy(const Vector& x)
{
    double fx = 0.0;
    double charge = 0;
    m_interface->updateGeometry(x);
    if (m_interface->HasNan()) {
        m_error = true;
        return 0;
    }
    fx = m_interface->CalculateEnergy(false);
    m_energy = fx;
    return fx;
}

double LBFGS::getEnergyGradient(const Vector& x)
{
    double fx = 0.0;
    double charge = 0;
    m_interface->updateGeometry(x);
    if (m_interface->HasNan()) {
        m_error = true;
        return 0;
    }
    fx = m_interface->CalculateEnergy(true);
    auto gradient = m_interface->Gradient();
    m_error = std::isnan(fx);

    for (int i = 0; i < m_atoms; ++i) {
        m_gradient[3 * i + 0] = gradient(i, 0) * (m_constrains[i]);
        m_gradient[3 * i + 1] = gradient(i, 1) * (m_constrains[i]);
        m_gradient[3 * i + 2] = gradient(i, 2) * (m_constrains[i]);
    }
    m_energy = fx;
    return fx;
}
