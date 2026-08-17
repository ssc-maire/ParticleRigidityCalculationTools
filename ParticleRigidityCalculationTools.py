import decimal as dec
import numpy as np
import pandas as pd

onekftinkm = 0.3048

protonRestMass = dec.Decimal(1.67262192e-27)  #kg
chargeOfElectron = dec.Decimal(1.60217663e-19) #C
c = dec.Decimal(299792458.0) #m/s

def allowForNonSeriesInputArgs(functionToModify):

    def functionWithGeneralInputArgs(*args, **kwargs):

        newArgs = []
        for inputArg in args:
            if type(inputArg) == pd.Series:
                newArg = inputArg
            elif type(inputArg) in [int, float]:
                newArg = pd.Series([inputArg])
            else:
                newArg = pd.Series(inputArg)
            newArgs.append(newArg)

        result = functionToModify(*newArgs, **kwargs)

        return result

    return functionWithGeneralInputArgs

def getAtomicMass(atomicNumber):

    A = [1.0,  4.0,  6.9,  9.0, 10.8, 12.0, 14.0, 16.0, 19.0, 20.2,\
       23.0, 24.3, 27.0, 28.1, 31.0, 32.1, 35.4, 39.9, 39.1, 40.1,\
       44.9, 47.9, 50.9, 52.0, 54.9, 55.8, 58.9, 58.7, 63.5, 65.4,\
       69.7, 72.6, 74.9, 79.0, 79.9, 83.8, 85.5, 87.6, 88.9, 91.2,\
       92.9, 95.9, 97.0,101.0,102.9,106.4,107.9,112.4,114.8,118.7,\
      121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1,140.9,144.2,\
      145.0,150.4,152.0,157.3,158.3,162.5,164.9,167.3,168.9,173.0,\
      175.0,178.5,180.9,183.9,186.2,190.2,192.2,195.1,197.0,200.6,\
      204.4,207.2,209.0,209.0,210.0,222.0,223.0,226.0,227.0,232.0,\
      231.0,238.0]

    if ( atomicNumber < 0 ):
        atomicMass = 1.0           # handle case for electron
    elif (atomicNumber <= 92):
        atomicMass = A[atomicNumber-1]       # look up table for other elements.
    else:
        atomicMass = 0.0           # just in case

    return atomicMass

def determineParticleAttributes(particleMassAU, particleChargeAU):
    m0 = dec.Decimal(particleMassAU) * protonRestMass #kg
    particleCharge = dec.Decimal(particleChargeAU) * chargeOfElectron #C
    particleRestEnergy = m0 * (c**2)
    return particleCharge,particleRestEnergy

def _as_series(values):
    if isinstance(values, pd.Series):
        return values
    if np.isscalar(values):
        return pd.Series([values])
    return pd.Series(values)

@allowForNonSeriesInputArgs
def convertParticleEnergyToRigidity(particleKineticEnergyInMeV:pd.Series, particleMassAU = 1, particleChargeAU = 1):
    """Convert total kinetic energy (MeV) to total rigidity (GV).

    ``particleKineticEnergyInMeV`` is the nucleus kinetic energy, not MeV/n.
    For a per-nucleon energy grid use ``convertPerNucleonEnergyToTotalRigidity``.
    """

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    particleKineticEnergyInJoules = particleKineticEnergyInMeV.apply(dec.Decimal) * chargeOfElectron * dec.Decimal(1e6)

    totalParticleEnergy = particleKineticEnergyInJoules + particleRestEnergy
    pc = np.sqrt((totalParticleEnergy**2) - (particleRestEnergy**2))

    #rigidity = pc / particleCharge

    rigidityInGV = (pc / dec.Decimal(particleCharge)) * dec.Decimal(1e-9)

    return rigidityInGV.apply(float)

@allowForNonSeriesInputArgs
def convertParticleRigidityToEnergy(particleRigidityInGV:pd.Series, particleMassAU = 1, particleChargeAU = 1):
    """Convert total rigidity (GV) to total kinetic energy (MeV).

    The returned energy is the nucleus kinetic energy, not MeV/n.
    For a per-nucleon energy grid use ``convertTotalRigidityToPerNucleonEnergy``.
    """

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    pc = particleRigidityInGV.apply(dec.Decimal) * particleCharge * dec.Decimal(1e9)

    totalParticleEnergy = np.sqrt((pc**2) + (particleRestEnergy**2))

    particleKEinJoules = totalParticleEnergy - particleRestEnergy

    KEinMeV = particleKEinJoules / (chargeOfElectron * dec.Decimal(1e6))

    return KEinMeV.apply(float)

def calculate_dKEoverdR(particleKineticEnergyInMeV:pd.Series, particleChargeInCoulombs, particleRestEnergy):
    particleKineticEnergyInJoules = particleKineticEnergyInMeV.apply(dec.Decimal) * chargeOfElectron * dec.Decimal(1e6)

    totalParticleEnergy = particleKineticEnergyInJoules + particleRestEnergy
    pc = np.sqrt((totalParticleEnergy**2) - (particleRestEnergy**2)) # units of pc are in joules

    #fullFactor = pc/particleKineticEnergyInJoules
    fullFactor = pc/totalParticleEnergy

    dKEInMeV_drigidityInGV = fullFactor * particleChargeInCoulombs * dec.Decimal(1e9) / (chargeOfElectron * dec.Decimal(1e6))
    return dKEInMeV_drigidityInGV # output units are in milli electron charges

@allowForNonSeriesInputArgs
def convertParticleEnergySpecToRigiditySpec(particleKineticEnergyInMeV:pd.Series, fluxInEnergyMeVform:pd.Series, particleMassAU = 1, particleChargeAU = 1):
    """Convert dJ/dE (per total MeV) to dJ/dR (per total GV).

    ``fluxInEnergyMeVform`` must be a differential flux per total MeV, not per MeV/n.
    For a per-nucleon spectrum use ``convertPerNucleonEnergySpecToTotalRigiditySpec``.
    """

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    dKEInMeV_drigidityInGV = calculate_dKEoverdR(particleKineticEnergyInMeV, particleCharge, particleRestEnergy)

    outputRigidities = convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = particleMassAU, particleChargeAU = particleChargeAU)
    outputRigiditySpectrum = (dKEInMeV_drigidityInGV * fluxInEnergyMeVform.apply(dec.Decimal)).apply(float)

    outputDataFrame = pd.DataFrame({"Rigidity":outputRigidities, "Rigidity distribution values":outputRigiditySpectrum})

    return outputDataFrame.map(float)

@allowForNonSeriesInputArgs
def convertParticleRigiditySpecToEnergySpec(particleRigidityInGV:pd.Series, fluxInRigidityGVform:pd.Series, particleMassAU = 1, particleChargeAU = 1):
    """Convert dJ/dR (per total GV) to dJ/dE (per total MeV).

    The returned energy column is total MeV and the flux is per total MeV.
    For a per-nucleon spectrum use ``convertTotalRigiditySpecToPerNucleonEnergySpec``.
    """

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    particleKineticEnergyInMeV = convertParticleRigidityToEnergy(particleRigidityInGV, particleMassAU = particleMassAU, particleChargeAU = particleChargeAU).apply(dec.Decimal)

    dKEInMeV_drigidityInGV = calculate_dKEoverdR(particleKineticEnergyInMeV, particleCharge, particleRestEnergy)

    outputEnergies = particleKineticEnergyInMeV

    dKEInMeV_drigidityInGV.replace(dec.Decimal(0),dec.Decimal(np.nan),inplace=True)
    outputEnergySpectrum = (fluxInRigidityGVform.apply(dec.Decimal) / dKEInMeV_drigidityInGV).apply(float)

    outputDataFrame = pd.DataFrame({"Energy":outputEnergies, "Energy distribution values":outputEnergySpectrum})

    return outputDataFrame.map(float)

def convertPerNucleonEnergyToTotalRigidity(
    energyPerNucleonMeV,
    particleMassAU=1,
    particleChargeAU=1,
):
    """Convert kinetic energy per nucleon (MeV/n) to total rigidity (GV).

    ``convertParticleEnergyToRigidity`` expects total kinetic energy in MeV, so
    this multiplies by mass number ``A`` first.
    """
    energy_per_nucleon = _as_series(energyPerNucleonMeV)
    energy_total_MeV = energy_per_nucleon * particleMassAU
    return convertParticleEnergyToRigidity(
        energy_total_MeV,
        particleMassAU=particleMassAU,
        particleChargeAU=particleChargeAU,
    )

def convertTotalRigidityToPerNucleonEnergy(
    particleRigidityInGV,
    particleMassAU=1,
    particleChargeAU=1,
):
    """Convert total rigidity (GV) to kinetic energy per nucleon (MeV/n)."""
    energy_total_MeV = convertParticleRigidityToEnergy(
        particleRigidityInGV,
        particleMassAU=particleMassAU,
        particleChargeAU=particleChargeAU,
    )
    return energy_total_MeV / particleMassAU

def convertPerNucleonEnergySpecToTotalRigiditySpec(
    energyPerNucleonMeV,
    fluxPerMeVPerNucleon,
    particleMassAU=1,
    particleChargeAU=1,
):
    """Convert a per-nucleon energy spectrum to a total-rigidity spectrum.

    If ``j_En`` is in (cm-2 s-1 sr-1 (MeV/n)-1) and ``R`` is total GV, then
    ``j_R = j_En * d(E/n)/dR = (j_En / A) * d(E_tot)/dR``.
    """
    energy_per_nucleon = _as_series(energyPerNucleonMeV)
    flux_per_mev_n = _as_series(fluxPerMeVPerNucleon)
    energy_total_MeV = energy_per_nucleon * particleMassAU
    flux_for_total_jacobian = flux_per_mev_n / particleMassAU
    return convertParticleEnergySpecToRigiditySpec(
        energy_total_MeV,
        flux_for_total_jacobian,
        particleMassAU=particleMassAU,
        particleChargeAU=particleChargeAU,
    )

def convertTotalRigiditySpecToPerNucleonEnergySpec(
    particleRigidityInGV,
    fluxInRigidityGVform,
    particleMassAU=1,
    particleChargeAU=1,
):
    """Convert a total-rigidity spectrum to a per-nucleon energy spectrum.

    If ``j_R`` is in (cm-2 s-1 sr-1 GV-1) and ``E_n`` is MeV/n, then
    ``j_En = j_R * dR/d(E/n) = A * j_Etot``.
    """
    total_energy_spec = convertParticleRigiditySpecToEnergySpec(
        particleRigidityInGV,
        fluxInRigidityGVform,
        particleMassAU=particleMassAU,
        particleChargeAU=particleChargeAU,
    )
    return pd.DataFrame({
        "Energy": total_energy_spec["Energy"] / particleMassAU,
        "Energy distribution values": total_energy_spec["Energy distribution values"] * particleMassAU,
    })

