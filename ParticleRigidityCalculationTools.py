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

@allowForNonSeriesInputArgs
def convertParticleEnergyToRigidity(particleKineticEnergyInMeV:pd.Series, particleMassAU = 1, particleChargeAU = 1):

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    particleKineticEnergyInJoules = particleKineticEnergyInMeV.apply(dec.Decimal) * chargeOfElectron * dec.Decimal(1e6)

    totalParticleEnergy = particleKineticEnergyInJoules + particleRestEnergy
    pc = np.sqrt((totalParticleEnergy**2) - (particleRestEnergy**2))

    #rigidity = pc / particleCharge

    rigidityInGV = (pc / dec.Decimal(particleCharge)) * dec.Decimal(1e-9)

    return rigidityInGV.apply(float)

@allowForNonSeriesInputArgs
def convertParticleRigidityToEnergy(particleRigidityInGV:pd.Series, particleMassAU = 1, particleChargeAU = 1):

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

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    dKEInMeV_drigidityInGV = calculate_dKEoverdR(particleKineticEnergyInMeV, particleCharge, particleRestEnergy)

    outputRigidities = convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = particleMassAU, particleChargeAU = particleChargeAU)
    outputRigiditySpectrum = (dKEInMeV_drigidityInGV * fluxInEnergyMeVform.apply(dec.Decimal)).apply(float)

    outputDataFrame = pd.DataFrame({"Rigidity":outputRigidities, "Rigidity distribution values":outputRigiditySpectrum})

    return outputDataFrame.applymap(float)

@allowForNonSeriesInputArgs
def convertParticleRigiditySpecToEnergySpec(particleRigidityInGV:pd.Series, fluxInRigidityGVform:pd.Series, particleMassAU = 1, particleChargeAU = 1):

    particleCharge, particleRestEnergy = determineParticleAttributes(particleMassAU, particleChargeAU)

    particleKineticEnergyInMeV = convertParticleRigidityToEnergy(particleRigidityInGV, particleMassAU = particleMassAU, particleChargeAU = particleChargeAU).apply(dec.Decimal)

    dKEInMeV_drigidityInGV = calculate_dKEoverdR(particleKineticEnergyInMeV, particleCharge, particleRestEnergy)

    outputEnergies = particleKineticEnergyInMeV

    dKEInMeV_drigidityInGV.replace(dec.Decimal(0),dec.Decimal(np.nan),inplace=True)
    outputEnergySpectrum = (fluxInRigidityGVform.apply(dec.Decimal) / dKEInMeV_drigidityInGV).apply(float)

    outputDataFrame = pd.DataFrame({"Energy":outputEnergies, "Energy distribution values":outputEnergySpectrum})

    return outputDataFrame.applymap(float)
