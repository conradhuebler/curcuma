/*
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/xtbinterface.h"

#include "src/capabilities/optimiser/LBFGSppInterface.h"

#include "src/capabilities/curcumaopt.h"

#include "src/core/molecule.h"

#include <string>

static const std::string coordinates = "C          -1.15550709190386    1.19752831253866    3.85081716354492\n"
                                       "C          -0.11105900681600    1.99614711803360    4.64249159516548\n"
                                       "C           1.26991467831331    1.38292683755582    4.40325206227649\n"
                                       "C           1.23404756229262   -0.10558940590441    4.78471249451125\n"
                                       "C           0.07011421404710   -0.80110219365967    4.06278873232617\n"
                                       "O          -1.16685528961610   -0.14234390233949    4.26757794868596\n"
                                       "O          -0.07504653131483    3.34345230406789    4.24289301857246\n"
                                       "O           2.24660986331273    2.00870313333009    5.19823074744100\n"
                                       "O           1.08904287305615   -0.24771856591227    6.17069133656958\n"
                                       "C          -0.12129952248181   -2.23672089376139    4.56705261848591\n"
                                       "O          -1.14590810084431   -2.89255902504769    3.87334668652737\n"
                                       "H          -0.98263239460084    3.67100326011242    4.24911610799025\n"
                                       "H           2.12546492416642    2.96122791063435    5.10004602734353\n"
                                       "H          -1.93683857475845   -2.34641750636149    3.95579035229608\n"
                                       "H          -0.36883832556162    1.90961680288030    5.70895134462519\n"
                                       "H           1.52073759215575    1.49006203621763    3.33604127135536\n"
                                       "H           2.17626644086572   -0.57409572175387    4.45548800920133\n"
                                       "H           0.29188729056438   -0.80941351557805    2.98156415202224\n"
                                       "H           1.68462987213601    0.39934647878840    6.57556583515199\n"
                                       "H           0.79198200862058   -2.80912635723135    4.39888657790903\n"
                                       "H          -0.32981627125567   -2.19114852097348    5.64196222880138\n"
                                       "O          -2.38679950876982    1.78743333965832    4.09018548846009\n"
                                       "H          -0.90806799422725    1.23936198777007    2.76409525795907\n"
                                       "C          -3.45612997442079    1.23193402738811    3.34809545989941\n"
                                       "H          -4.33088577244393    1.83023419069101    3.59173900208989\n"
                                       "H          -3.25100629990380    1.28939166422417    2.27332005412277\n"
                                       "H          -3.62736380667261    0.19103332835231    3.62892675772879\n";

int main()
{
    int GFNmethod = 1;
    /* Empty molecule */
    Molecule molecule;

    /* just load the xyz coordindates given as simple string */
    molecule.LoadMolecule(coordinates);

    /* print what just got loaded */
    molecule.print_geom();

    /* Using CurcumaOpt::LBFGSOptimise function, which gives some console output and so on */
    /* It is controlled by a json structure */
    static json CurcumaOptJson{
        { "Solver", 1 },
        { "writeXYZ", true },
        { "printOutput", true },
        { "dE", 0.1 },
        { "dRMSD", 0.01 },
        { "GFN", GFNmethod },
        { "InnerLoop", 20 },
        { "OuterLoop", 100 },
        { "LBFGS_eps", 1e-5 },
        { "Threads", 1 },
        { "Charge", 0 },
        { "Spin", 0 }
    };
    std::string output;
    std::vector<Molecule> intermediate;
    Molecule final = CurcumaOpt::LBFGSOptimise(&molecule, CurcumaOptJson, output, &intermediate);
    final.print_geom();

    /* That was it for dealing with the function */
    /* Using the CurcumaOpt class, several molecules can be optimised, in parallel even */
    /* For that have a look at src/main.cpp after
     * strcmp(argv[1], "-opt") == 0 */

    /* The interface can be called directly as follows */

    XTBInterface interface;
    interface.InitialiseMolecule(molecule);

    //double final_energy = interface.GFNCalculation(1); // 2 = GFN2

    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 20; // or any number of LBFGS iterations

    LBFGSSolver<double> solver(param); // LBFGS solver

    /* This interface converts the molecule to xtb readable structures and calls energy and gradient calculation */

    LBFGSInterface fun(3 * molecule.AtomCount());
    fun.setMolecule(&molecule);
    fun.setInterface(&interface);
    fun.setMethod(GFNmethod);

    Vector parameter(3 * molecule.AtomCount());
    Geometry geometry = molecule.getGeometry(); // It is just an Eigen::MatrixXd
    for (int i = 0; i < molecule.AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }
    double fx;
    int niter = solver.minimize(fun, parameter, fx); // No more control over optimisation, just wait till end

    parameter = fun.Parameter(); // Get optimised coordinates

    for (int i = 0; i < molecule.AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    Molecule result2 = molecule;
    result2.setGeometry(geometry);
    result2.setEnergy(fun.LastEnergy()); // read last energy

    result2.print_geom();

    return 0;
}
