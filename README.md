# ParticleRigidityCalculationTools

A set of Python functions for processing and converting data in units of particle kinetic energy to units of particle rigidity (typically expressed in GV; gigavolts) and vice versa. This is frequently done in the field of solar system radiation physics.

# Installation

To install from pypi, run

```
sudo pip3 install ParticleRigidityCalculationTools
```
Alternatively, you can install directly from the Github repository.

To do this you can first clone the repository, and then from the cloned repository, run

```
pip install .
```

(Use `sudo pip install .` for a system-wide install.)

# Usage

For all the functions contained in this module, **total** kinetic energy is expressed in **MeV** (megaelectronvolts), and **total** rigidity is expressed in **GV** (gigavolts) unless otherwise stated.

The low-level converters `convertParticleEnergyToRigidity` and `convertParticleRigidityToEnergy` take the nucleus kinetic energy, **not** MeV/n. Passing a per-nucleon energy grid for Z>1 without multiplying by mass number \(A\) underestimates rigidity (for helium, typically by a factor of about two). Use `convertPerNucleonEnergyToTotalRigidity` / `convertTotalRigidityToPerNucleonEnergy` when the energy grid is in MeV/n, which is the usual convention for heavy-ion spectra and for MAIRE-style atmospheric response matrices.

## General Conversion Functions

To convert particle kinetic energy to rigidity, use the `convertParticleEnergyToRigidity` function. 

Particle kinetic energies can be supplied as a float, int, list, [NumPy array](https://numpy.org/doc/stable/reference/generated/numpy.array.html) or [Pandas Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html). Particle mass in atomic units should also be supplied, as well as the particle charge magnitude in atomic units, as particle rigidity is dependent on these quantities.

For instance, to calculate particle rigidities for several kinetic energies at once you can first define a list of particle kinetic energies:
```
import ParticleRigidityCalculationTools as PRCT

particleKineticEnergyInMeV = [250.0, 578.5, 1056.8, 5123.9]
```

and then running
```
PRCT.convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = 1.0, particleChargeAU = 1.0)
```

will give the corresponding rigidities for a **proton** with kinetic energies of 250.0 MeV, 578.5 MeV, 1056.8 MeV and 5123.9 MeV respectively:
```
0    0.729134
1    1.191740
2    1.760670
3    5.989121
```
note that the output to this function, as with all rigidity calculation functions in this module is a [Pandas Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html). 

To perform the opposite calculation, calculating kinetic energies from a list of rigidities, you can use the `convertParticleRigidityToEnergy` function, which uses exactly the same input format but using input rigidities instead of energies. Using the output from the previous function:

```
outputtedRigiditiesSeries = PRCT.convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = 1.0, particleChargeAU = 1.0)
```
we can get back the original set of proton kinetic energies with
```
PRCT.convertParticleRigidityToEnergy(outputtedRigiditiesSeries,particleMassAU=1.0,particleChargeAU=1.0)
```
which returns
```
0     250.0
1     578.5
2    1056.8
3    5123.9
```
as a [Pandas Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html).

When not using protons, you can either directly input the particle mass from tabulated values or use the `getAtomicMass` function to get tabulated average mass values for a particle with a particular atomic number. For instance, for an alpha particle/helium ion:

```
alphaParticleAtomicNumber = 2

PRCT.getAtomicMass(alphaParticleAtomicNumber)
```

returns
```
4.0
```

The charge magnitude `|Z|` used by all functions is the particle atomic number (for electrons, `|Z| = 1`). `particleChargeAU` may be signed; only the magnitude is used, because magnetic rigidity is \(R = pc/|q|\).

For a helium ion at 1129 MeV/n, the total rigidity is obtained from the **total** kinetic energy \(A \times 1129\) MeV:

```
PRCT.convertPerNucleonEnergyToTotalRigidity(1129.0, particleMassAU=4.0, particleChargeAU=2)
```

which returns about **3.69 GV**. Calling `convertParticleEnergyToRigidity(1129.0, particleMassAU=4.0, particleChargeAU=2)` treats 1129 as total MeV and incorrectly returns about 1.56 GV.

## Electrons

Electrons are not nuclei. Their rest energy is \(mc^{2} \approx 0.511\,\mathrm{MeV}\), not the proton rest energy \(\approx 938\,\mathrm{MeV}\). Pass the electron-to-proton mass ratio as `particleMassAU`, and charge magnitude `|Z| = 1`.

`getAtomicMass(-1)` returns that mass ratio \(m_e/m_p \approx 1/1836\), using the CODATA 2018 electron and proton masses. Do **not** use `particleMassAU = 1`: that treats the electron as a proton and overestimates rigidity at MeV energies by a factor of about 30.

```
electronMassAU = PRCT.getAtomicMass(-1)

PRCT.convertParticleEnergyToRigidity(1.0, particleMassAU=electronMassAU, particleChargeAU=1)
```

returns about **0.00142 GV** (1.42 MV) for a 1 MeV electron. That matches the relativistic definition

\[
R = \frac{\sqrt{K(K+2mc^{2})}}{|Z|\times 1000}\ \mathrm{GV}
\]

with \(K\) and \(mc^{2}\) in MeV. `particleChargeAU=-1` gives the same positive rigidity, because only \(|q|\) is used.

At high energy the electron is ultrarelativistic, so \(R \approx (K + 0.511)/1000\) GV. At 10 GeV that is about **10.0005 GV**.

Use the total-energy converters (`convertParticleEnergyToRigidity`, `convertParticleEnergySpecToRigiditySpec`, and their inverses) with the electron mass. Do **not** use the per-nucleon wrappers for electrons: those multiply energy by mass number \(A\), which is the right correction for ions such as helium, not for \(m_e/m_p\).

## Spectrum Conversion Functions

A user might not necessarily want to just convert individual numbers between units of rigidity and energy, they might also want to convert a kinetic energy distribution or rigidity distribution. This might usually be expressed in the form of $\frac{dN}{dE}$ or $\frac{dN}{dR}$, where E and R are particle kinetic energy and rigidity respectively, and where both quantities are expressed in terms of kinetic energy and rigidity respectively. As there is a one-to-one relationship between kinetic energy and rigidity, it is possible to analytically convert between these two quantities using $\frac{dN}{dR} = \frac{dN}{dE} \times \frac{dE}{dR}$, where $\frac{dR}{dE}$ can be calculated using the definition of the [magnetic rigidity of a particle](https://www.nmdb.eu/public_outreach/de/07_md/).

Tools are available in this module to perform all of this process automatically. The function `convertParticleEnergySpecToRigiditySpec` can be used to convert kinetic energy distributions into rigidity distributions, for example:

```
energyValuesInMeV = [1000,2000,3000,4000,5000]
energyDistributionValues = [1,0.5,0.2,0.1,0.01]

PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 1.0,particleChargeAU = 1.0)
```

returns
```
   Rigidity  Rigidity distribution values
0  1.696038                    875.025647
1  2.784437                    473.822152
2  3.824870                    194.241037
3  4.848317                     98.178407
4  5.863678                      9.874384
```
as a [Pandas DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html).

The function `convertParticleRigiditySpecToEnergySpec` can be used to perform the opposite operation, converting particle rigidity to kinetic energy. For example, 

```
rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 1,particleChargeAU = 1)

PRCT.convertParticleRigiditySpecToEnergySpec(rigiditySpec["Rigidity"],rigiditySpec["Rigidity distribution values"],particleMassAU = 1,particleChargeAU = 1)
```

returns
```
   Energy  Energy distribution values
0  1000.0                        1.00
1  2000.0                        0.50
2  3000.0                        0.20
3  4000.0                        0.10
4  5000.0                        0.01
```
the original kinetic energies and distribution values that were used for the energy distribution.

The same round-trip for a **per-nucleon** helium spectrum uses `convertPerNucleonEnergySpecToTotalRigiditySpec` and `convertTotalRigiditySpecToPerNucleonEnergySpec`. If \(j_{E_n}\) is in particles cm\(^{-2}\) s\(^{-1}\) sr\(^{-1}\) (MeV/n)\(^{-1}\), the rigidity flux is \(j_R = j_{E_n}\,\mathrm{d}(E/n)/\mathrm{d}R\).
