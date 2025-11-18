/*
 * Test compilation of CV framework
 * Claude Generated - November 2025
 */

#include "src/capabilities/cv/collective_variable.h"
#include "src/capabilities/cv/cv_distance.h"
#include "src/capabilities/cv/cv_angle.h"
#include "src/capabilities/cv/cv_dihedral.h"
#include "src/capabilities/cv/cv_gyration.h"
#include "src/capabilities/cv/cv_coordination.h"
#include "src/capabilities/cv/cv_factory.h"
#include "src/capabilities/cv/bias_engine.h"

int main() {
    // Test CV creation via factory
    auto dist_cv = CV::CVFactory::create("distance", {0, 10});
    auto angle_cv = CV::CVFactory::create("angle", {5, 0, 10});
    auto dih_cv = CV::CVFactory::create("dihedral", {1, 2, 3, 4});
    auto rg_cv = CV::CVFactory::create("gyration", {});
    auto coord_cv = CV::CVFactory::create("coordination", {0}, {1, 5, 9});

    // Test BiasEngine
    CV::BiasEngine engine;
    engine.addCV(CV::CVFactory::create("distance", {0, 10}), 0.2);
    engine.setWellTempered(true, 3000.0);
    engine.setInitialHeight(0.1);

    return 0;
}
