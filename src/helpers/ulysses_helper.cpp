#include "Eigen/Dense"

#include "src/core/interface/ulysses.h"

#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Geometry;

int main()
{
    std::vector<double> geom = {
        1.49574345397488, 1.16281825573451, 0.28849847338774,
        -0.49945625222174, -2.53040659067913, -0.53071797893969,
        0.16194648100926, -0.87962369198586, 0.35253076320581,
        3.22412736186075, -2.24117494192853, -0.20154564030214,
        3.97051978668812, -1.07886850522647, -0.23108198001996,
        3.34784608607901, 0.14765154772948, -0.06232336424791,
        1.97897334251832, 0.20833216808745, 0.14161901458194,
        1.22521319238394, -0.95170242664151, 0.16076903226233,
        1.84402933789503, -2.18703464982174, -0.01791521590997,
        1.09952904099103, -3.47995606086188, 0.00132142847995,
        -0.21243863808497, -3.42313523883092, -0.29781606419764,
        1.62179782894737, -4.52960069465670, 0.25829369761659,
        3.92945856088550, 1.05722760598657, -0.08135305033683,
        5.03838404064105, -1.12730466913327, -0.38177929265508,
        3.68742240287454, -3.20834081773612, -0.32281097212030
    };
    int num_atoms = 15;
    Geometry matrix = Eigen::MatrixXd::Zero(num_atoms, 3);

    for (int i = 0; i < num_atoms; ++i) {
        matrix(i, 0) = geom[i * 3];
        matrix(i, 1) = geom[i * 3 + 1];
        matrix(i, 2) = geom[i * 3 + 2];
    }

    std::vector<int> atoms = { 1, 1, 1, 6, 6, 6, 6, 6, 6, 6, 8, 8, 1, 1, 1 };
    int charge = 0;
    int spin = 0;
    int mult = 1;
    {
        std::string method = "gfn2";
        UlyssesObject ulysses;
        ulysses.setMethod(method);

        ulysses.setMolecule(matrix, atoms, charge, spin, "C1");
        ulysses.Calculate(true, true);
        std::cout << ulysses.Energy() << std::endl;
        std::cout << ulysses.Gradient() << std::endl;
    }

    {
        std::string method = "pm6-d3h4x";
        UlyssesObject ulysses;
        ulysses.setMethod(method);

        ulysses.setMolecule(matrix, atoms, charge, mult, "C1");
        ulysses.Calculate(true, true);
        std::cout << ulysses.Energy() << std::endl;
        std::cout << ulysses.Gradient() << std::endl;
    }
    return 0;
}