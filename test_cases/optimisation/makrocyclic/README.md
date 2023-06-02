# Overview
The structure files AnGrad.xyz and NumGrad.xyz were obtained using the git commit 4e19c68f8c252babdcd7c87e95287c8f34177e9a of curcuma and commit 7e3848617795ddd0e25f4b772e679adfee583229 of the LBFGSpp submodule.
The best result (lowest energy) was obtained using analytical gradients, although curcuma ended with
```sh
LBFGS interface signalled some logic error!
 --         the moving direction increases the objective function value         --
```
. Using the numerical gradients, the energy remains higher, but no error is signaled. I can not claim to have a perfect implementation of the universal force field, but at least MD simulation, which need a gradient work quite nice now.

Using the latest LBFGSpp changes (master branch of my fork), the energy is higher than with the older linesearch implementations. Numerical gradients are even better than analytical.

## Old implementation

|NumGrad: |185          |3.027176        |0.000000       | 0.000000     |   0.028869     |   0.000000
|AnGrad : |365       |   3.015478     |   0.000000   |     0.067949       | 0.015598     |   0.061000

## New implementation

Can easily be tested using the AppImage (https://github.com/conradhuebler/curcuma/releases/download/0.0.108/curcuma-nightly-0.0.108-20_04-x86_64-Linux.AppImage)
For analytical gradients
```sh
curcuma-nightly-0.0.108-20_04-x86_64-Linux.AppImage -opt input.xyz -gradient 0
```
or for numerical gradients
```sh
curcuma-nightly-0.0.108-20_04-x86_64-Linux.AppImage -opt input.xyz -gradient 1
```
|NumGrad: |147     |     3.027298   |     -1.063203    |   0.150695   |     0.031905    |    0.311000
|AnGrad : |138     |     3.028958   |    -24.520050     |  0.149910    |    0.046107    |    0.063000

So, analytical gradients with the old implementation signal some error, but result in the lowest energy (3.015478), every else result is higher. Similar stuff happen with other combination of potential methods (besides UFF eq. GFNFF, GFN2), but having a look at UFF would be the easiest.
It would be fine to catch some exceptions etc within curcuma and signal some optimisations as success. But the new implementation results in a too high last energy.
