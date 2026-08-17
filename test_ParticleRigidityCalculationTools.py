
import numpy as np
import pandas as pd
import ParticleRigidityCalculationTools as PRCT

def test_getAtomicMass():

    assert PRCT.getAtomicMass(12) == 24.3

def test_rigidityConversion():

    particleKineticEnergyInMeV = [250.0, 578.5, 1056.8, 5123.9]

    outputValues = PRCT.convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = 1.0, particleChargeAU = 1.0)

    roundedOutputValues = outputValues.apply(lambda x:round(x,2))

    assert list(roundedOutputValues) == [0.73, 1.19,1.76, 5.99]

def oneDPround(inputVal):

    return round(inputVal,1)

def test_rigidityAndEnergyConversion():

    particleKineticEnergyInMeV = [250.0, 578.5, 1056.8, 5123.9]

    outputRigidityValues = PRCT.convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = 1.0, particleChargeAU = 1.0)

    outputEnergyValues = PRCT.convertParticleRigidityToEnergy(outputRigidityValues, particleMassAU = 1.0, particleChargeAU = 1.0)

    assert list(outputEnergyValues.apply(oneDPround)) == particleKineticEnergyInMeV


#######################################################

def test_EnergySpecToRigiditySpec():

    energyValuesInMeV = [1000, 2000, 3000, 4000, 5000]
    energyDistributionValues = [1, 0.5, 0.2, 0.1, 0.01]

    rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 1,particleChargeAU = 1)

    assert rigiditySpec["Rigidity"].round(2).iloc[2] == 3.82
    assert rigiditySpec["Rigidity distribution values"].round(2).iloc[2] == 194.24

def test_RigiditySpecToEnergySpec():

    rigidityValuesInGV = np.linspace(0.0,10.0,19)
    rigidityDistributionValues = np.linspace(0.0,10.0,19)

    energySpec = PRCT.convertParticleRigiditySpecToEnergySpec(rigidityValuesInGV,rigidityDistributionValues,particleMassAU = 1,particleChargeAU = 1)

def test_BothSpectralConversions():

    energyValuesInMeV = [1000, 2000, 3000, 4000, 5000]
    energyDistributionValues = [1, 0.5, 0.2, 0.1, 0.01]

    rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 1,particleChargeAU = 1)

    energySpec = PRCT.convertParticleRigiditySpecToEnergySpec(rigiditySpec["Rigidity"],rigiditySpec["Rigidity distribution values"],particleMassAU = 1,particleChargeAU = 1)

    assert energySpec["Energy"].round(2).array == pd.Series(energyValuesInMeV).apply(float).round(2).array
    assert energySpec["Energy distribution values"].round(2).array == pd.Series(energyDistributionValues).apply(float).round(2).array


def test_BothSpectralConversionsAlpha():

    energyValuesInMeV = [1000, 2000, 3000, 4000, 5000]
    energyDistributionValues = [1, 0.5, 0.2, 0.1, 0.01]

    rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 4,particleChargeAU = 2)

    energySpec = PRCT.convertParticleRigiditySpecToEnergySpec(rigiditySpec["Rigidity"],rigiditySpec["Rigidity distribution values"],particleMassAU = 4,particleChargeAU = 2)

    assert energySpec["Energy"].round(2).array == pd.Series(energyValuesInMeV).apply(float).round(2).array
    assert energySpec["Energy distribution values"].round(2).array == pd.Series(energyDistributionValues).apply(float).round(2).array


def test_helium_per_nucleon_energy_is_not_total_energy():
    energy_per_nucleon_MeV = 1129.0
    correct_rigidity = PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_per_nucleon_MeV,
        particleMassAU=4.0,
        particleChargeAU=2,
    )
    incorrect_rigidity = PRCT.convertParticleEnergyToRigidity(
        energy_per_nucleon_MeV,
        particleMassAU=4.0,
        particleChargeAU=2,
    )
    np.testing.assert_allclose(correct_rigidity.iloc[0], 3.684, rtol=1e-3)
    np.testing.assert_allclose(incorrect_rigidity.iloc[0], 1.561, rtol=1e-3)
    assert correct_rigidity.iloc[0] > 2.0 * incorrect_rigidity.iloc[0]


def test_helium_per_nucleon_rigidity_is_twice_proton_at_same_energy_per_nucleon():
    energy_per_nucleon_MeV = [10.0, 100.0, 1000.0, 1129.0]
    proton_rigidity = PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_per_nucleon_MeV,
        particleMassAU=1.0,
        particleChargeAU=1,
    )
    helium_rigidity = PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_per_nucleon_MeV,
        particleMassAU=4.0,
        particleChargeAU=2,
    )
    np.testing.assert_allclose(helium_rigidity, 2.0 * proton_rigidity, rtol=1e-12)


def test_per_nucleon_energy_and_rigidity_round_trip():
    energy_per_nucleon_MeV = [250.0, 578.5, 1056.8, 1129.0]
    rigidity = PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_per_nucleon_MeV,
        particleMassAU=4.0,
        particleChargeAU=2,
    )
    recovered_energy = PRCT.convertTotalRigidityToPerNucleonEnergy(
        rigidity,
        particleMassAU=4.0,
        particleChargeAU=2,
    )
    np.testing.assert_allclose(recovered_energy, energy_per_nucleon_MeV, rtol=1e-10)


def test_per_nucleon_spectral_conversions_round_trip():
    energy_per_nucleon_MeV = [1000, 2000, 3000, 4000, 5000]
    flux_per_mev_n = [1, 0.5, 0.2, 0.1, 0.01]

    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        energy_per_nucleon_MeV,
        flux_per_mev_n,
        particleMassAU=4,
        particleChargeAU=2,
    )
    energy_spec = PRCT.convertTotalRigiditySpecToPerNucleonEnergySpec(
        rigidity_spec["Rigidity"],
        rigidity_spec["Rigidity distribution values"],
        particleMassAU=4,
        particleChargeAU=2,
    )

    np.testing.assert_allclose(energy_spec["Energy"], energy_per_nucleon_MeV, rtol=1e-10)
    np.testing.assert_allclose(energy_spec["Energy distribution values"], flux_per_mev_n, rtol=1e-10)

    naive_rigidity_spec = PRCT.convertParticleEnergySpecToRigiditySpec(
        energy_per_nucleon_MeV,
        flux_per_mev_n,
        particleMassAU=4,
        particleChargeAU=2,
    )
    assert not np.allclose(
        rigidity_spec["Rigidity"],
        naive_rigidity_spec["Rigidity"],
        rtol=1e-2,
    )
    np.testing.assert_allclose(
        rigidity_spec["Rigidity"],
        2.0 * PRCT.convertParticleEnergyToRigidity(energy_per_nucleon_MeV, particleMassAU=1, particleChargeAU=1),
        rtol=1e-12,
    )
