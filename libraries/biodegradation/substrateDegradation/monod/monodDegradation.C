/*---------------------------------------------------------------------------*\
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "monodDegradation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace landfillBiodegradationModels
    {
        defineTypeNameAndDebug(monodDegradation, 0);

        addToRunTimeSelectionTable
        (
            landfillBiodegradation,
            monodDegradation,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::landfillBiodegradationModels::monodDegradation::monodDegradation
(
    const word& name,
    const dictionary& biodegradProperties,
    volScalarField& substrate
)
    :
    landfillBiodegradation(name, biodegradProperties, substrate),
    monodCoeffs_(biodegradProperties.subDict(typeName + "Coeffs")),
    biomassA0_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("biomassA0", dimensionedScalar("biomassA0", dimDensity, 0.15))),
    biomassN0_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("biomassN0", dimensionedScalar("biomassN0", dimDensity, 0.15))),
    biomassA_
    (
        IOobject
        (
            "biomassA",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        substrate.mesh(),
        biomassA0_
    ),
    biomassN_
    (
        IOobject
        (
            "biomassN",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        substrate.mesh(),
        biomassN0_
    ),
    inertMass_
    (
        IOobject
        (
            "inertMass",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        substrate.mesh(),
        dimensionedScalar(dimDensity,0.0)
    ),
    molConcO2_(substrate_.db().lookupObjectRef<volScalarField>("molConcentrationO2")),
    T_(substrate_.db().lookupObject<volScalarField>("T")),
    p_(substrate_.db().lookupObject<volScalarField>("p")),
    stoicO2A_(monodCoeffs_.lookupOrDefault<scalar>("stoicO2A", -6)),
    stoicO2N_(monodCoeffs_.lookupOrDefault<scalar>("stoicO2N", 0)),
    stoicCO2A_(monodCoeffs_.lookupOrDefault<scalar>("stoicCO2A", 6)),
    stoicCO2N_(monodCoeffs_.lookupOrDefault<scalar>("stoicCO2N", 3)),
    stoicCH4A_(monodCoeffs_.lookupOrDefault<scalar>("stoicCH4A", 0)),
    stoicCH4N_(monodCoeffs_.lookupOrDefault<scalar>("stoicCH4N", 3)),
    stoicH2OA_(monodCoeffs_.lookupOrDefault<scalar>("stoicH2OA", 5)),
    stoicH2ON_(monodCoeffs_.lookupOrDefault<scalar>("stoicH2ON", -1)),
    kMaxA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("kMaxA", dimensionedScalar("kMaxA", dimless/dimTime, 1.157e-5))),
    kMaxN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("kMaxN", dimensionedScalar("kMaxN", dimless/dimTime, 2.315e-7))),
    KSA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("KSA", dimensionedScalar("KSA", dimDensity, 181))),
    KSN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("KSN", dimensionedScalar("KSN", dimDensity, 46))),
    KO2_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("KO2", dimensionedScalar("KO2", dimMoles/dimVolume, 0.8173))),
    YSA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("YSA", dimensionedScalar("YSA", dimless, 0.3))),
    YSN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("YSN", dimensionedScalar("YSN", dimless, 0.05))),
    productFactorA_((1-YSA_)/YSA_),
    productFactorN_((1-YSN_)/YSN_),
    molWeightSubstrate_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("molWeightSubstrate", dimensionedScalar("molWeightSubstrate", dimMass/dimMoles, 0.162))), // kg/mol
    molConcO2Anaerobic_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("molConcO2Anaerobic", dimensionedScalar("molConcO2Anaerobic", dimMoles/dimVolume, 0.040))), // mol/m3
    decayConstant_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("decayConstant", dimensionedScalar("decayConstant", dimless, 0.05))),
    enthalpyA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("enthalpyA", dimensionedScalar("enthalpyA", dimEnergy/dimMass, 1.736e7))), // J/kg-substrate
    enthalpyN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("enthalpyN", dimensionedScalar("enthalpyN", dimEnergy/dimMass, 1.672e6))), // J/kg-substrate
    TminA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("TminA", dimensionedScalar("TminA", dimTemperature, 278.15))),
    TmaxA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("TmaxA", dimensionedScalar("TmaxA", dimTemperature, 344.75))),
    ToptA_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("ToptA", dimensionedScalar("ToptA", dimTemperature, 331.75))),
    TminN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("TminN", dimensionedScalar("TminN", dimTemperature, 288.15))),
    TmaxN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("TmaxN", dimensionedScalar("TmaxN", dimTemperature, 338.15))),
    ToptN_(monodCoeffs_.lookupOrDefault<dimensionedScalar>("ToptN", dimensionedScalar("ToptN", dimTemperature, 318.15))),
    kTempA_
    (
        IOobject
        (
            "kTempA",
            T_.time().timeName(),
            T_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar(dimless,0.0)
    ),
    kTempN_
    (
        IOobject
        (
            "kTempN",
            T_.time().timeName(),
            T_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar(dimless,0.0)
    ),
    growthA_
    (
        IOobject
        (
            "growthA",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        substrate_.mesh(),
        dimensionedScalar(dimDensity/dimTime,0.0)
    ),
    growthN_
    (
        IOobject
        (
            "growthN",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        substrate_.mesh(),
        dimensionedScalar(dimDensity/dimTime,0.0)
    ),
    decayA_
    (
        IOobject
        (
            "decayA",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        substrate_.mesh(),
        dimensionedScalar(dimDensity/dimTime,0.0)
    ),
    decayN_
    (
        IOobject
        (
            "decayN",
            substrate_.time().timeName(),
            substrate_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        substrate_.mesh(),
        dimensionedScalar(dimDensity/dimTime,0.0)
    )
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * *//
void Foam::landfillBiodegradationModels::monodDegradation::updatekTempA_()
{
    /*forAll(T_,celli)
    {
        if (T_[celli] < TminA_.value() || T_[celli] > TmaxA_.value())
        {
            kTempA_[celli] = 0;
        }
        else
        {
            kTempA_[celli] = (T_[celli] - TmaxA_.value())*pow((T_[celli] - TminA_.value()),2)
                /(ToptA_.value() - TminA_.value())
                /(((ToptA_.value() - TminA_.value())*(T_[celli] - ToptA_.value())-(ToptA_.value() - TmaxA_.value())*(ToptA_.value()+TminA_.value()-2*T_[celli]))+VSMALL);
        }
    }*/

    //write with pos()
    kTempA_ = pos(T_ - TminA_) * pos(TmaxA_ - T_) * (T_ - TmaxA_) * pow(T_ - TminA_,2)
        /(ToptA_ - TminA_)
        /(((ToptA_ - TminA_)*(T_ - ToptA_)-(ToptA_ - TmaxA_)*(ToptA_+TminA_-2*T_))+dimensionedScalar(dimTemperature*dimTemperature,VSMALL));
}

void Foam::landfillBiodegradationModels::monodDegradation::updatekTempN_()
{
    /*forAll(T_,celli)
    {
        if (T_[celli] <= TminN_.value() || T_[celli] >= TmaxN_.value())
        {
            kTempN_[celli] = 0;
        }
        else
        {
            kTempN_[celli] = (T_[celli] - TmaxN_.value())*pow((T_[celli] - TminN_.value()),2)
                /(ToptN_.value() - TminN_.value())
                /(((ToptN_.value() - TminN_.value())*(T_[celli] - ToptN_.value())-(ToptN_.value() - TmaxN_.value())*(ToptN_.value()+TminN_.value()-2*T_[celli]))+VSMALL);
        }
    }*/

    //write with pos()
    kTempN_ = pos(T_ - TminN_) * pos(TmaxN_ - T_) * (T_ - TmaxN_) * pow(T_ - TminN_,2)
        /(ToptN_ - TminN_)
        /(((ToptN_ - TminN_)*(T_ - ToptN_)-(ToptN_ - TmaxN_)*(ToptN_+TminN_-2*T_))+dimensionedScalar(dimTemperature*dimTemperature,VSMALL));
}

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::kSat_(volScalarField& s, dimensionedScalar& k)
{
    tmp<volScalarField> tkSat
    (
        new volScalarField
        (
            IOobject
            (
                "kSat",
                s.time().timeName(),
                s.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            s.mesh(),
            dimensionedScalar(dimless,0.0)
        )
    );

    volScalarField& kSat = tkSat.ref();

    kSat = pos(s) * s/(s+k);

    return tkSat;
}

void Foam::landfillBiodegradationModels::monodDegradation::updategrowthA()
{
    growthA_ = kMaxA_ * kTempA_ * kSat_(substrate_, KSA_) * kSat_(molConcO2_,KO2_) * biomassA_;

}

void Foam::landfillBiodegradationModels::monodDegradation::updategrowthN()
{

    growthN_ = pos(molConcO2Anaerobic_ - molConcO2_) * kMaxN_ * kTempN_ * kSat_(substrate_, KSN_) * biomassN_;


}

void Foam::landfillBiodegradationModels::monodDegradation::updatedecayA()
{
    decayA_ = - decayConstant_ * kMaxA_  * (biomassA_ - biomassA0_);
}

void Foam::landfillBiodegradationModels::monodDegradation::updatedecayN()
{
    decayN_ = - decayConstant_ * kMaxN_  * (biomassN_ - biomassN0_);
}

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourceSubstrate()
{
    tmp<volScalarField> sourceSubstrate
    (
        new volScalarField
        (
            IOobject
            (
                "sourceSubstrate",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimDensity/dimTime,0.0)
        )
    );

    sourceSubstrate = - growthA_ / YSA_
        - growthN_ / YSN_;


    return sourceSubstrate;
}

/*Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourceSubstrateAerobic()
{
    tmp<volScalarField> sourceSubstrateAerobic
    (
        new volScalarField
        (
            IOobject
            (
                "sourceSubstrateAerobic",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimDensity/dimTime,0.0)
        )
    );

    sourceSubstrateAerobic = - growthA() / YSA_;

    return sourceSubstrateAerobic;
}*/

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourceO2()
{
    tmp<volScalarField> sourceO2
    (
        new volScalarField
        (
            IOobject
            (
                "sourceO2",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimMoles/dimVolume/dimTime,0.0)
        )
    );

    sourceO2 = stoicO2A_  / molWeightSubstrate_ * productFactorA_ * growthA_
        + stoicO2N_  / molWeightSubstrate_ * productFactorN_ * growthN_;

    return sourceO2;
}

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourceCO2()
{
    tmp<volScalarField> sourceCO2
    (
        new volScalarField
        (
            IOobject
            (
                "sourceCO2",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimMoles/dimVolume/dimTime,0.0)
        )
    );

    sourceCO2 = stoicCO2A_  / molWeightSubstrate_ * productFactorA_ * growthA_
        + stoicCO2N_  / molWeightSubstrate_ * productFactorN_ * growthN_;

    return sourceCO2;
}

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourceCH4()
{
    tmp<volScalarField> sourceCH4
    (
        new volScalarField
        (
            IOobject
            (
                "sourceCH4",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimMoles/dimVolume/dimTime,0.0)
        )
    );


    sourceCH4 = stoicCH4A_  / molWeightSubstrate_ * productFactorA_ * growthA_
        + stoicCH4N_  / molWeightSubstrate_ * productFactorN_ * growthN_;

    return sourceCH4;
}

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourcew()
{
    tmp<volScalarField> sourcew
    (
        new volScalarField
        (
            IOobject
            (
                "sourcew",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimDensity/dimTime,0.0)
        )
    );

    dimensionedScalar molwH2O("molwH2O", dimMass/dimMoles, 0.018);

    sourcew = stoicH2OA_ * molwH2O / molWeightSubstrate_ * productFactorA_ * growthA_
        + stoicH2ON_ * molwH2O / molWeightSubstrate_ * productFactorN_ * growthN_;

    return sourcew;
}

Foam::tmp<Foam::volScalarField> Foam::landfillBiodegradationModels::monodDegradation::sourceHeat()
{
    tmp<volScalarField> sourceHeat
    (
        new volScalarField
        (
            IOobject
            (
                "sourceHeat",
                substrate_.time().timeName(),
                substrate_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            substrate_.mesh(),
            dimensionedScalar(dimEnergy/dimTime,0.0)
        )
    );

    sourceHeat = enthalpyA_ * growthA_ * productFactorA_
        + enthalpyN_ * growthN_ * productFactorN_;

    return sourceHeat;
}

void Foam::landfillBiodegradationModels::monodDegradation::update()
{
    updatekTempA_();
    updatekTempN_();
    updategrowthA();
    updategrowthN();
    updatedecayA();
    updatedecayN();

    fvScalarMatrix biomassAEqn
    (
        fvm::ddt(biomassA_)
        ==
        decayA_
        + growthA_
    );
    fvScalarMatrix biomassNEqn
    (
        fvm::ddt(biomassN_)
        ==
        decayN_
        + growthN_
    );
    fvScalarMatrix substrateEqn
    (
        fvm::ddt(substrate_)
        ==
        sourceSubstrate()
    );
    fvScalarMatrix inertMassEqn
    (
        fvm::ddt(inertMass_)
        ==
        -decayA_ - decayN_
    );

    biomassAEqn.solve();
    biomassNEqn.solve();
    substrateEqn.solve();
    inertMassEqn.solve();

    Info << "Maximum biomassA = " << gMax(biomassA_.internalField()) << endl;
    Info << "Minimum biomassA = " << gMin(biomassA_.internalField()) << endl;
    Info << "Maximum biomassN = " << gMax(biomassN_.internalField()) << endl;
    Info << "Minimum biomassN = " << gMin(biomassN_.internalField()) << endl;
    Info << "Maximum inertMass = " << gMax(inertMass_.internalField()) << endl;
    Info << "Minimum inertMass = " << gMin(inertMass_.internalField()) << endl;
}

