# ColorRCMC
Solution of Radiative Transfer Equation in a silica medium with spherical plain and core-shell nanoparticles.

Works on GPU with MATLAB

Execute main_core_shell_core_shell_mix.m to run the code for core-shell core-shell particle mix.
Execute main_core_shell_core_shell_mix.m to run the code for core-shell plain particle mix.
To consider single type of particle rather than a mixture, set volume fraction of 2nd type to 0.

The code first calculates size dependent dielectric functions. Then analytical Mie solution for a single core-shell and/or plain nanoparticle is performed to estimate effective medium properties. Afterwards, it calls monte carlo method to estimate reflectance, finally net radiative cooling and color of the coating are estimated at post process.

You can contact me from refetali@gmail.com

In case of use, please cite the article as:, Colored Radiative Cooling Coatings with Nanoparticles, ACS Photonics, Refet Ali Yalçın, Etienne Blandre, Karl Joulain, Jeremie Drevillon, https://doi.org/10.1021/acsphotonics.0c00513
