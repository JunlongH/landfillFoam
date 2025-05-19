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

#include "methaneOxidation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::methaneOxidation::methaneOxidation
(
        const dictionary& reactionProperties,
        const volScalarField& T,
        const volScalarField& p,
        const volScalarField& Sw
)
:
    methaneOxidationCoeffs_(reactionProperties.subDict("methaneOxidation")),
    T_(T),
    p_(p),
    Sw_(Sw),
    molConcCH4_(T_.db().lookupObject<volScalarField>("molConcentrationCH4")),
    molConcO2_(T_.db().lookupObject<volScalarField>("molConcentrationO2")),
    stoicCH4_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicCH4", -1.0)),
    stoicO2_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicO2", -2.0)),
    stoicCO2_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicCO2", 1.0)),
    stoicH2O_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicH2O", 2.0)),
    Swmin_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("Swmin",0.3)),
    Swmax_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("Swmax",0.99)),
    Swfc_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("Swfc",0.5)),
    kmax_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("kmax", dimensionedScalar("kmax",dimDensity/dimTime,2.4e-6))),
    KmCH4_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("KmCH4", dimensionedScalar("KmCH4",dimMoles/dimVolume,0.2697))),
    KmO2_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("KmO2", dimensionedScalar("KmO2",dimMoles/dimVolume,0.4904))),
    enthalpyCH4_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("enthalpyCH4",dimensionedScalar("enthalpyCH4",dimEnergy/dimMass,6.956e6))),
    kTemp_
    (
        IOobject
        (
            "kTemp",
            T_.time().timeName(),
            T_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar(dimless,1.0)
    ),
    kSw_
    (
        IOobject
        (
            "kSw",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        dimensionedScalar(dimless,1.0)
    ),
    kO2_
    (
        IOobject
        (
            "kO2",
            molConcO2_.time().timeName(),
            molConcO2_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        molConcO2_.mesh(),
        dimensionedScalar(dimless,1.0)
    ),
    kCH4_
    (
        IOobject
        (
            "kCH4",
            molConcCH4_.time().timeName(),
            molConcCH4_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        molConcCH4_.mesh(),
        dimensionedScalar(dimless,1.0)
    )
{}

Foam::methaneOxidation::methaneOxidation
(
        const dictionary& reactionProperties,
        const volScalarField& T,
        const volScalarField& p,
        const volScalarField& Sw,
        const volScalarField& molConcCH4,
        const volScalarField& molConcO2
)
:
    methaneOxidationCoeffs_(reactionProperties.subDict("methaneOxidation")),
    T_(T),
    p_(p),
    Sw_(Sw),
    molConcCH4_(molConcCH4),
    molConcO2_(molConcO2),
    stoicCH4_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicCH4", -1.0)),
    stoicO2_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicO2", -2.0)),
    stoicCO2_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicCO2", 1.0)),
    stoicH2O_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("stoicH2O", 2.0)),
    Swmin_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("Swmin",0.3)),
    Swmax_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("Swmax",0.99)),
    Swfc_(methaneOxidationCoeffs_.lookupOrDefault<scalar>("Swfc",0.5)),
    kmax_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("kmax", dimensionedScalar("kmax",dimDensity/dimTime,2.4e-6))),
    KmCH4_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("KmCH4", dimensionedScalar("KmCH4",dimMoles/dimVolume,0.2697))),
    KmO2_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("KmO2", dimensionedScalar("KmO2",dimMoles/dimVolume,0.4904))),
    enthalpyCH4_(methaneOxidationCoeffs_.lookupOrDefault<dimensionedScalar>("enthalpyCH4",dimensionedScalar("enthalpyCH4",dimEnergy/dimMass,6.956e6))),
    kTemp_
    (
        IOobject
        (
            "kTemp",
            T_.time().timeName(),
            T_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar(dimless,1.0)
    ),
    kSw_
    (
        IOobject
        (
            "kSw",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        dimensionedScalar(dimless,1.0)
    ),
    kO2_
    (
        IOobject
        (
            "kO2",
            molConcO2_.time().timeName(),
            molConcO2_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        molConcO2_.mesh(),
        dimensionedScalar(dimless,1.0)
    ),
    kCH4_
    (
        IOobject
        (
            "kCH4",
            molConcCH4_.time().timeName(),
            molConcCH4_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        molConcCH4_.mesh(),
        dimensionedScalar(dimless,1.0)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * *//
void Foam::methaneOxidation::updatekTemp_()
{
    
    // alternate scalar temperature into dimensionedScalar
    kTemp_ = pos(T_ - dimensionedScalar(dimTemperature,273.15)) * pos(dimensionedScalar(dimTemperature,288.15) - T_) * dimensionedScalar(dimless/dimTemperature, 0.0142) * (T_ - dimensionedScalar(dimTemperature,273.15))
        + pos(T_ - dimensionedScalar(dimTemperature,288.15)) * pos(dimensionedScalar(dimTemperature,306.15) - T_) * (dimensionedScalar(dimless/dimTemperature, 0.112)  * (T_ - dimensionedScalar(dimTemperature,273.15)) - 1.47)
        + pos(T_ - dimensionedScalar(dimTemperature,306.15)) * pos(dimensionedScalar(dimTemperature,318.15) - T_) * (2.235 - (dimensionedScalar(dimless/dimTemperature, 0.18) * (T_ - dimensionedScalar(dimTemperature,306.15))));

    /*forAll(T_,celli)
    {
        if(T_[celli] > 273.15 && T_[celli] < 288.15)
        {
            kTemp_[celli] = 0.0142 * (T_[celli] - 273.15);
        }
        else if(T_[celli] >= 288.15 && T_[celli] < 306.15)
        {
            kTemp_[celli] = 0.112 * (T_[celli] - 273.15) - 1.47;
        }
        else if(T_[celli] >= 306.15 && T_[celli] < 318.15)
        {
            kTemp_[celli] = 2.235 - 0.18 * (T_[celli] - 306.15);
        }
        else
        {
            kTemp_[celli] = 0;
        }
    }*/
}

void Foam::methaneOxidation::updatekSw_()
{
    kSw_ = pos(Sw_ - Swmin_) * pos(Swfc_ - Sw_) * (Sw_ - Swmin_) / (Swfc_ - Swmin_)
        + pos(Sw_ - Swfc_) * pos(Swmax_ - Sw_) * 1.0;

    /*forAll(Sw_,celli)
    {
        if(Sw_[celli] > Swmin_ && Sw_[celli] < Swfc_)
        {
            kSw_[celli] = (Sw_[celli] - Swmin_)/(Swfc_ - Swmin_);
        }
        else if(Sw_[celli] >= Swfc_ && Sw_[celli] < Swmax_)
        {
            kSw_[celli] = 1.0;
        }
        else
        {
            kSw_[celli] = 0;
        }
    }*/
}

void Foam::methaneOxidation::updatekO2_()
{
    kO2_ = pos(molConcO2_) * molConcO2_ / (molConcO2_ + KmO2_);
    /*forAll(kO2_,celli)
    {
        if(molConcO2_[celli] > 0.0)
        {
            kO2_[celli] = molConcO2_[celli]/(molConcO2_[celli] + KmO2_.value());
        }
        else
        {
            kO2_[celli] = 0.0;
        }
    }*/
}

void Foam::methaneOxidation::updatekCH4_()
{
    kCH4_ = pos(molConcCH4_) * molConcCH4_ / (molConcCH4_ + KmCH4_);
    /*forAll(kCH4_,celli)
    {
        if(molConcCH4_[celli] > 0.0)
        {
            kCH4_[celli] = molConcCH4_[celli]/(molConcCH4_[celli] + KmCH4_.value());
        }
        else
        {
            kCH4_[celli] = 0.0;
        }
    }*/
}

Foam::tmp<Foam::volScalarField> Foam::methaneOxidation::sourceCH4()
{
    tmp<volScalarField> tsourceCH4
    (
        new volScalarField
        (
            IOobject
            (
                "sourceCH4",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimMoles/dimVolume/dimTime
        )
    );

    dimensionedScalar molwCH4("molwCH4",dimMass/dimMoles,0.01604);

    tsourceCH4 = stoicCH4_ * kmax_ * kTemp_ * kSw_ * kO2_ * kCH4_ / molwCH4;

    return tsourceCH4;
}

Foam::tmp<Foam::volScalarField> Foam::methaneOxidation::sourceO2()
{
    tmp<volScalarField> tsourceO2
    (
        new volScalarField
        (
            IOobject
            (
                "sourceO2",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimMoles/dimVolume/dimTime
        )
    );

    dimensionedScalar molwCH4("molwCH4",dimMass/dimMoles,0.01604);

    tsourceO2 = stoicO2_ * kmax_ * kTemp_ * kSw_ * kO2_ * kCH4_ / molwCH4;

    return tsourceO2;
}

Foam::tmp<Foam::volScalarField> Foam::methaneOxidation::sourceCO2()
{
    tmp<volScalarField> tsourceCO2
    (
        new volScalarField
        (
            IOobject
            (
                "sourceCO2",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimMoles/dimVolume/dimTime
        )
    );

    dimensionedScalar molwCH4("molwCH4",dimMass/dimMoles,0.01604);

    tsourceCO2 = stoicCO2_ * kmax_ * kTemp_ * kSw_ * kO2_ * kCH4_ / molwCH4;

    return tsourceCO2;
}

Foam::tmp<Foam::volScalarField> Foam::methaneOxidation::sourcew()
{
    tmp<volScalarField> tsourcew
    (
        new volScalarField
        (
            IOobject
            (
                "sourcew",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimMoles/dimVolume/dimTime
        )
    );

    dimensionedScalar molwCH4("molwCH4",dimMass/dimMoles,0.01604);
    dimensionedScalar molwH2O("molwH2O",dimMass/dimMoles,0.018015);

    tsourcew = stoicH2O_ * kmax_ * kTemp_ * kSw_ * kO2_ * kCH4_ / molwCH4 * molwH2O;

    return tsourcew;
}

Foam::tmp<Foam::volScalarField> Foam::methaneOxidation::sourceHeat()
{
    tmp<volScalarField> tsourceHeat
    (
        new volScalarField
        (
            IOobject
            (
                "sourceHeat",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimEnergy/dimVolume/dimTime
        )
    );

    tsourceHeat = enthalpyCH4_ * kmax_ * kTemp_ * kSw_ * kO2_ * kCH4_;

    return tsourceHeat;
}

void Foam::methaneOxidation::update()
{
    updatekTemp_();
    updatekSw_();
    updatekO2_();
    updatekCH4_();
}


