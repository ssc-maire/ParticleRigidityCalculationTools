
import numpy as np
import ParticleRigidityCalculationTools as PRCT

# Tabulated atomic masses for proton, helium, magnesium, and electrons (m_e / m_p).
def test_getAtomicMass():

    assert PRCT.getAtomicMass(1) == 1.0
    assert PRCT.getAtomicMass(2) == 4.0
    assert PRCT.getAtomicMass(12) == 24.3
    np.testing.assert_allclose(PRCT.getAtomicMass(-1), ELECTRON_MASS_AU, rtol=1e-12)

# Proton total kinetic energy (MeV) to total rigidity (GV) against known values.
def test_rigidityConversion():

    particleKineticEnergyInMeV = [250.0, 578.5, 1056.8, 5123.9]

    outputValues = PRCT.convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = 1.0, particleChargeAU = 1.0)

    np.testing.assert_allclose(
        outputValues,
        [0.7291337629108413, 1.1917395085827414, 1.7606697947406509, 5.989121464609327],
        rtol=1e-12,
    )

# Proton energy -> rigidity -> energy returns the original total kinetic energies.
def test_rigidityAndEnergyConversion():

    particleKineticEnergyInMeV = [250.0, 578.5, 1056.8, 5123.9]

    outputRigidityValues = PRCT.convertParticleEnergyToRigidity(particleKineticEnergyInMeV, particleMassAU = 1.0, particleChargeAU = 1.0)

    outputEnergyValues = PRCT.convertParticleRigidityToEnergy(outputRigidityValues, particleMassAU = 1.0, particleChargeAU = 1.0)

    np.testing.assert_allclose(outputEnergyValues, particleKineticEnergyInMeV, rtol=1e-10)


#######################################################

# Proton dJ/dE (per total MeV) to dJ/dR (per total GV) against known values.
def test_EnergySpecToRigiditySpec():

    energyValuesInMeV = [1000, 2000, 3000, 4000, 5000]
    energyDistributionValues = [1, 0.5, 0.2, 0.1, 0.01]

    rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 1,particleChargeAU = 1)

    np.testing.assert_allclose(
        rigiditySpec["Rigidity"],
        [1.6960377875702215, 2.78443681087077, 3.8248702632374703, 4.848316894290674, 5.863678102038896],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        rigiditySpec["Rigidity distribution values"],
        [875.0256466528116, 473.8221524534995, 194.2410365434809, 98.1784074969649, 9.874384357464102],
        rtol=1e-12,
    )

# Proton dJ/dR (per total GV) to dJ/dE (per total MeV); R=0 maps to zero energy and NaN flux.
def test_RigiditySpecToEnergySpec():

    rigidityValuesInGV = np.linspace(0.0,10.0,19)
    rigidityDistributionValues = np.linspace(0.0,10.0,19)

    energySpec = PRCT.convertParticleRigiditySpecToEnergySpec(rigidityValuesInGV,rigidityDistributionValues,particleMassAU = 1,particleChargeAU = 1)

    expected_energy = PRCT.convertParticleRigidityToEnergy(
        rigidityValuesInGV, particleMassAU=1, particleChargeAU=1
    )
    np.testing.assert_allclose(energySpec["Energy"], expected_energy, rtol=1e-10)
    np.testing.assert_allclose(energySpec["Energy"].iloc[0], 0.0, atol=1e-15)
    assert np.isnan(energySpec["Energy distribution values"].iloc[0])
    np.testing.assert_allclose(energySpec["Energy"].iloc[9], 4149.001691493931, rtol=1e-12)
    np.testing.assert_allclose(
        energySpec["Energy distribution values"].iloc[9],
        0.005087273779926977,
        rtol=1e-12,
    )

    recovered_rigidity_spec = PRCT.convertParticleEnergySpecToRigiditySpec(
        energySpec["Energy"].iloc[1:],
        energySpec["Energy distribution values"].iloc[1:],
        particleMassAU=1,
        particleChargeAU=1,
    )
    np.testing.assert_allclose(recovered_rigidity_spec["Rigidity"], rigidityValuesInGV[1:], rtol=1e-10)
    np.testing.assert_allclose(
        recovered_rigidity_spec["Rigidity distribution values"],
        rigidityDistributionValues[1:],
        rtol=1e-10,
    )

# Proton total-energy and total-rigidity spectra convert into each other and back.
def test_BothSpectralConversions():

    energyValuesInMeV = [1000, 2000, 3000, 4000, 5000]
    energyDistributionValues = [1, 0.5, 0.2, 0.1, 0.01]

    rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 1,particleChargeAU = 1)

    energySpec = PRCT.convertParticleRigiditySpecToEnergySpec(rigiditySpec["Rigidity"],rigiditySpec["Rigidity distribution values"],particleMassAU = 1,particleChargeAU = 1)

    np.testing.assert_allclose(energySpec["Energy"], energyValuesInMeV, rtol=1e-10)
    np.testing.assert_allclose(energySpec["Energy distribution values"], energyDistributionValues, rtol=1e-10)


# Helium total-MeV (not MeV/n) energy and rigidity spectra convert into each other and back.
def test_BothSpectralConversionsAlpha():

    energyValuesInMeV = [1000, 2000, 3000, 4000, 5000]
    energyDistributionValues = [1, 0.5, 0.2, 0.1, 0.01]

    rigiditySpec = PRCT.convertParticleEnergySpecToRigiditySpec(energyValuesInMeV,energyDistributionValues,particleMassAU = 4,particleChargeAU = 2)

    energySpec = PRCT.convertParticleRigiditySpecToEnergySpec(rigiditySpec["Rigidity"],rigiditySpec["Rigidity distribution values"],particleMassAU = 4,particleChargeAU = 2)

    np.testing.assert_allclose(energySpec["Energy"], energyValuesInMeV, rtol=1e-10)
    np.testing.assert_allclose(energySpec["Energy distribution values"], energyDistributionValues, rtol=1e-10)


# Helium 1129 MeV/n must use the per-nucleon wrapper (~3.68 GV), not total-MeV conversion (~1.56 GV).
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
    np.testing.assert_allclose(correct_rigidity.iloc[0], 3.6841603524720896, rtol=1e-12)
    np.testing.assert_allclose(incorrect_rigidity.iloc[0], 1.561178601468076, rtol=1e-12)
    assert correct_rigidity.iloc[0] > 2.0 * incorrect_rigidity.iloc[0]


# At the same energy per nucleon, helium total rigidity is exactly twice the proton rigidity (A/Z = 2).
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


# Helium MeV/n -> total GV -> MeV/n recovers the original per-nucleon energies.
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


# Helium per-nucleon spectrum round-trip, and the naive total-MeV converter must disagree on both R and dJ/dR.
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
    assert not np.allclose(
        rigidity_spec["Rigidity distribution values"],
        naive_rigidity_spec["Rigidity distribution values"],
        rtol=1e-2,
    )
    np.testing.assert_allclose(
        rigidity_spec["Rigidity"],
        2.0 * PRCT.convertParticleEnergyToRigidity(energy_per_nucleon_MeV, particleMassAU=1, particleChargeAU=1),
        rtol=1e-12,
    )


HELIUM_A = 4.0
HELIUM_Z = 2
PROTON_A = 1.0
PROTON_Z = 1

# Same physical constants as ParticleRigidityCalculationTools, used only for an
# independent float64 check of the Decimal conversion chain.
_PROTON_REST_MASS_KG = 1.67262192e-27
_ELECTRON_REST_MASS_KG = 9.1093837015e-31
_ELECTRON_CHARGE_C = 1.60217663e-19
_SPEED_OF_LIGHT_M_S = 299792458.0
ELECTRON_MASS_AU = _ELECTRON_REST_MASS_KG / _PROTON_REST_MASS_KG
ELECTRON_CHARGE_MAGNITUDE_AU = 1
CODATA_ELECTRON_REST_ENERGY_MEV = 0.51099895

HELIUM_PER_NUCLEON_ENERGY_MEV = [11.294627058970837, 28.370820458389794, 1129.0]
HELIUM_PER_NUCLEON_RIGIDITY_GV_KGO = [
    0.2920440736449699,
    0.4649473147684023,
    3.6841603524720896,
]
HELIUM_PER_NUCLEON_FLUX = [5.565730138613413e-6, 1.3394663729029165e-5, 1.0]
HELIUM_PER_NUCLEON_RIGIDITY_FLUX_KGO = [
    4.279421540294369e-4,
    1.6106808615034475e-3,
    445.5340413443853,
]
HELIUM_SPECTRUM_ENERGY_MEV_N = [1000, 2000, 3000, 4000, 5000]
HELIUM_SPECTRUM_FLUX_PER_MEV_N = [1, 0.5, 0.2, 0.1, 0.01]
HELIUM_SPECTRUM_RIGIDITY_GV_KGO = [
    3.392075575140443,
    5.56887362174154,
    7.6497405264749405,
    9.696633788581348,
    11.727356204077791,
]
HELIUM_SPECTRUM_RIGIDITY_FLUX_KGO = [
    437.5128233264058,
    236.91107622674974,
    97.12051827174045,
    49.08920374848245,
    4.937192178732051,
]


# Independent float64 R = pc/|q| from total kinetic energy, used to check the Decimal converters.
def _independent_total_rigidity_gv(ke_total_mev, mass_au, charge_magnitude_au):
    rest_energy_j = mass_au * _PROTON_REST_MASS_KG * (_SPEED_OF_LIGHT_M_S ** 2)
    kinetic_energy_j = np.asarray(ke_total_mev, dtype=np.float64) * _ELECTRON_CHARGE_C * 1e6
    total_energy_j = kinetic_energy_j + rest_energy_j
    pc_j = np.sqrt(total_energy_j ** 2 - rest_energy_j ** 2)
    return (pc_j / (abs(charge_magnitude_au) * _ELECTRON_CHARGE_C)) * 1e-9


# Independent float64 dE_tot/dR (MeV/GV), used to check the helium dJ/dR Jacobian.
def _independent_d_total_energy_d_rigidity(ke_total_mev, mass_au, charge_magnitude_au):
    rest_energy_j = mass_au * _PROTON_REST_MASS_KG * (_SPEED_OF_LIGHT_M_S ** 2)
    charge_magnitude_c = abs(charge_magnitude_au) * _ELECTRON_CHARGE_C
    kinetic_energy_j = np.asarray(ke_total_mev, dtype=np.float64) * _ELECTRON_CHARGE_C * 1e6
    total_energy_j = kinetic_energy_j + rest_energy_j
    pc_j = np.sqrt(total_energy_j ** 2 - rest_energy_j ** 2)
    return (pc_j / total_energy_j) * charge_magnitude_c * 1e9 / (_ELECTRON_CHARGE_C * 1e6)


# Helium MeV/n -> total GV against high-precision known values and the independent analytic formula.
def test_helium_per_nucleon_rigidity_matches_kgo():
    rigidity = PRCT.convertPerNucleonEnergyToTotalRigidity(
        HELIUM_PER_NUCLEON_ENERGY_MEV,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    )
    independent_rigidity = _independent_total_rigidity_gv(
        np.asarray(HELIUM_PER_NUCLEON_ENERGY_MEV) * HELIUM_A,
        HELIUM_A,
        HELIUM_Z,
    )

    np.testing.assert_allclose(rigidity, HELIUM_PER_NUCLEON_RIGIDITY_GV_KGO, rtol=1e-12)
    np.testing.assert_allclose(rigidity, independent_rigidity, rtol=1e-12)


# Helium dJ/d(E/n) -> dJ/dR against known values and the independent analytic Jacobian.
def test_helium_per_nucleon_rigidity_spectrum_matches_kgo():
    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        HELIUM_PER_NUCLEON_ENERGY_MEV,
        HELIUM_PER_NUCLEON_FLUX,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    )
    independent_jacobian = _independent_d_total_energy_d_rigidity(
        np.asarray(HELIUM_PER_NUCLEON_ENERGY_MEV) * HELIUM_A,
        HELIUM_A,
        HELIUM_Z,
    ) / HELIUM_A
    independent_rigidity_flux = np.asarray(HELIUM_PER_NUCLEON_FLUX) * independent_jacobian

    np.testing.assert_allclose(
        rigidity_spec["Rigidity"],
        HELIUM_PER_NUCLEON_RIGIDITY_GV_KGO,
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        rigidity_spec["Rigidity distribution values"],
        HELIUM_PER_NUCLEON_RIGIDITY_FLUX_KGO,
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        rigidity_spec["Rigidity distribution values"],
        independent_rigidity_flux,
        rtol=1e-12,
    )


# Helium per-nucleon spectrum on the 1000–5000 MeV/n grid against high-precision R and dJ/dR values.
def test_helium_per_nucleon_spectrum_grid_matches_kgo():
    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        HELIUM_SPECTRUM_ENERGY_MEV_N,
        HELIUM_SPECTRUM_FLUX_PER_MEV_N,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    )

    np.testing.assert_allclose(
        rigidity_spec["Rigidity"],
        HELIUM_SPECTRUM_RIGIDITY_GV_KGO,
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        rigidity_spec["Rigidity distribution values"],
        HELIUM_SPECTRUM_RIGIDITY_FLUX_KGO,
        rtol=1e-12,
    )


# At the same E/n and dJ/d(E/n), helium dJ/dR is half the proton dJ/dR because dR/d(E/n) is twice as large.
def test_helium_rigidity_flux_is_half_proton_flux_at_same_energy_per_nucleon():
    energy_per_nucleon_MeV = [10.0, 100.0, 1000.0, 1129.0]
    flux_per_mev_n = [1.0, 0.5, 0.2, 0.1]

    proton_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        energy_per_nucleon_MeV,
        flux_per_mev_n,
        particleMassAU=PROTON_A,
        particleChargeAU=PROTON_Z,
    )
    helium_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        energy_per_nucleon_MeV,
        flux_per_mev_n,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    )

    np.testing.assert_allclose(
        helium_spec["Rigidity distribution values"],
        0.5 * proton_spec["Rigidity distribution values"],
        rtol=1e-12,
    )


# Helium unit-flux dJ/dR equals the finite-difference d(E/n)/dR at 1129 MeV/n.
def test_helium_per_nucleon_jacobian_matches_finite_difference():
    energy_per_nucleon_MeV = 1129.0
    delta_energy = 1e-4
    rigidity_plus = PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_per_nucleon_MeV + delta_energy,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    ).iloc[0]
    rigidity_minus = PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_per_nucleon_MeV - delta_energy,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    ).iloc[0]
    d_energy_per_nucleon_d_rigidity = (2.0 * delta_energy) / (rigidity_plus - rigidity_minus)

    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        [energy_per_nucleon_MeV],
        [1.0],
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    )

    np.testing.assert_allclose(
        rigidity_spec["Rigidity distribution values"].iloc[0],
        d_energy_per_nucleon_d_rigidity,
        rtol=1e-8,
    )


def _electron_rest_energy_mev():
    return _ELECTRON_REST_MASS_KG * (_SPEED_OF_LIGHT_M_S ** 2) / (_ELECTRON_CHARGE_C * 1e6)


def _textbook_rigidity_gv(kinetic_energy_mev, rest_energy_mev, charge_magnitude_au):
    """R = pc/|q| = sqrt(K(K+2mc^2)) / (|Z| * 1000) GV, with K and mc^2 in MeV."""
    kinetic_energy_mev = np.asarray(kinetic_energy_mev, dtype=np.float64)
    pc_mev = np.sqrt(kinetic_energy_mev * (kinetic_energy_mev + 2.0 * rest_energy_mev))
    return pc_mev / (abs(charge_magnitude_au) * 1e3)


def _textbook_dKE_dR_mev_per_gv(kinetic_energy_mev, rest_energy_mev, charge_magnitude_au):
    """dK/dR = |Z| * 1000 * (pc) / E, with E = K + mc^2 in MeV and R in GV."""
    kinetic_energy_mev = np.asarray(kinetic_energy_mev, dtype=np.float64)
    total_energy_mev = kinetic_energy_mev + rest_energy_mev
    pc_mev = np.sqrt(kinetic_energy_mev * (kinetic_energy_mev + 2.0 * rest_energy_mev))
    return (abs(charge_magnitude_au) * 1e3) * pc_mev / total_energy_mev


# Electron rest energy from PRCT's mass and c matches CODATA mc^2 ≈ 0.511 MeV.
def test_electron_rest_energy_matches_codata():
    np.testing.assert_allclose(
        _electron_rest_energy_mev(),
        CODATA_ELECTRON_REST_ENERGY_MEV,
        rtol=1e-8,
    )


# Electron rigidity uses R = sqrt(K(K+2mc^2)) / 1000 GV with mc^2 = 0.511 MeV, not proton mass.
def test_electron_energy_to_rigidity_matches_textbook_formula():
    kinetic_energy_mev = [0.1, 0.511, 1.0, 10.0, 100.0, 1000.0]
    electron_mass_au = PRCT.getAtomicMass(-1)
    rigidity = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=ELECTRON_CHARGE_MAGNITUDE_AU,
    )
    textbook_rigidity = _textbook_rigidity_gv(
        kinetic_energy_mev,
        _electron_rest_energy_mev(),
        ELECTRON_CHARGE_MAGNITUDE_AU,
    )
    independent_si_rigidity = _independent_total_rigidity_gv(
        kinetic_energy_mev,
        electron_mass_au,
        ELECTRON_CHARGE_MAGNITUDE_AU,
    )

    np.testing.assert_allclose(rigidity, textbook_rigidity, rtol=1e-12)
    np.testing.assert_allclose(rigidity, independent_si_rigidity, rtol=1e-12)
    np.testing.assert_allclose(rigidity.iloc[2], 0.0014219697263106033, rtol=1e-12)

    proton_rigidity = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=PROTON_A,
        particleChargeAU=PROTON_Z,
    )
    assert rigidity.iloc[2] < 0.05 * proton_rigidity.iloc[2]


# Ultrarelativistic electrons satisfy R ≈ (K + mc^2) / 1000 GV.
def test_electron_ultrarelativistic_rigidity_limit():
    kinetic_energy_mev = 10000.0
    rigidity = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=PRCT.getAtomicMass(-1),
        particleChargeAU=ELECTRON_CHARGE_MAGNITUDE_AU,
    ).iloc[0]
    np.testing.assert_allclose(
        rigidity,
        (kinetic_energy_mev + _electron_rest_energy_mev()) / 1000.0,
        rtol=1e-8,
    )


# Electron energy <-> rigidity round-trip with m_e/m_p and |Z| = 1.
def test_electron_energy_rigidity_round_trip():
    kinetic_energy_mev = [0.1, 1.0, 10.0, 100.0, 1000.0]
    electron_mass_au = PRCT.getAtomicMass(-1)
    rigidity = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=ELECTRON_CHARGE_MAGNITUDE_AU,
    )
    recovered_energy = PRCT.convertParticleRigidityToEnergy(
        rigidity,
        particleMassAU=electron_mass_au,
        particleChargeAU=ELECTRON_CHARGE_MAGNITUDE_AU,
    )
    np.testing.assert_allclose(recovered_energy, kinetic_energy_mev, rtol=1e-10)


# Magnetic rigidity uses charge magnitude |q|, so a negative particleChargeAU matches the positive value.
def test_charge_sign_is_ignored():
    kinetic_energy_mev = [1.0, 100.0]
    electron_mass_au = PRCT.getAtomicMass(-1)
    positive = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=1,
    )
    negative = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=-1,
    )
    np.testing.assert_allclose(negative, positive, rtol=1e-12)
    assert (positive > 0).all()

    helium_positive = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=HELIUM_A,
        particleChargeAU=HELIUM_Z,
    )
    helium_negative = PRCT.convertParticleEnergyToRigidity(
        kinetic_energy_mev,
        particleMassAU=HELIUM_A,
        particleChargeAU=-HELIUM_Z,
    )
    np.testing.assert_allclose(helium_negative, helium_positive, rtol=1e-12)

    flux_per_mev = [1.0, 0.5]
    positive_spec = PRCT.convertParticleEnergySpecToRigiditySpec(
        kinetic_energy_mev,
        flux_per_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=1,
    )
    negative_spec = PRCT.convertParticleEnergySpecToRigiditySpec(
        kinetic_energy_mev,
        flux_per_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=-1,
    )
    np.testing.assert_allclose(
        negative_spec["Rigidity"],
        positive_spec["Rigidity"],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        negative_spec["Rigidity distribution values"],
        positive_spec["Rigidity distribution values"],
        rtol=1e-12,
    )


# Electron dJ/dE -> dJ/dR uses j_R = j_E * |Z| * 1000 * (pc) / E.
def test_electron_energy_spectrum_jacobian_matches_textbook_formula():
    kinetic_energy_mev = [0.1, 1.0, 10.0, 100.0, 1000.0]
    flux_per_mev = [1.0, 0.5, 0.2, 0.1, 0.01]
    electron_mass_au = PRCT.getAtomicMass(-1)
    rigidity_spec = PRCT.convertParticleEnergySpecToRigiditySpec(
        kinetic_energy_mev,
        flux_per_mev,
        particleMassAU=electron_mass_au,
        particleChargeAU=ELECTRON_CHARGE_MAGNITUDE_AU,
    )
    textbook_flux = np.asarray(flux_per_mev) * _textbook_dKE_dR_mev_per_gv(
        kinetic_energy_mev,
        _electron_rest_energy_mev(),
        ELECTRON_CHARGE_MAGNITUDE_AU,
    )

    np.testing.assert_allclose(
        rigidity_spec["Rigidity"],
        _textbook_rigidity_gv(
            kinetic_energy_mev,
            _electron_rest_energy_mev(),
            ELECTRON_CHARGE_MAGNITUDE_AU,
        ),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        rigidity_spec["Rigidity distribution values"],
        textbook_flux,
        rtol=1e-12,
    )

    energy_spec = PRCT.convertParticleRigiditySpecToEnergySpec(
        rigidity_spec["Rigidity"],
        rigidity_spec["Rigidity distribution values"],
        particleMassAU=electron_mass_au,
        particleChargeAU=ELECTRON_CHARGE_MAGNITUDE_AU,
    )
    np.testing.assert_allclose(energy_spec["Energy"], kinetic_energy_mev, rtol=1e-10)
    np.testing.assert_allclose(
        energy_spec["Energy distribution values"],
        flux_per_mev,
        rtol=1e-10,
    )
