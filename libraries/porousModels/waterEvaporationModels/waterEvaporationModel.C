/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of landfill modelling project, an extension of OpenFOAM
    modeling landfill aeration.

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
#include "waterEvaporationModel.H"
#include "constants.H"
#include "binaryDiffusionCoeff.H"
#include "objectRegistry.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::waterEvaporationModel::waterEvaporationModel
(
    const word evapModel,
    const dimensionedScalar kEvap,
    const dimensionedScalar rho,
    const dimensionedScalar Vm,
    const dimensionedScalar as,
    const dimensionedScalar lDiff,
    const volScalarField& Sw,
    const volScalarField& eps,
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& pc,
    const volScalarField& volFracw
)
    :
    evapModel_(evapModel),
    kEvap_(kEvap),
    rho_(rho),
    Vm_(Vm),
    as_(as),
    lDiff_(lDiff),
    Sw_(Sw),
    eps_(eps),
    T_(T),
    p_(p),
    pc_(pc),
    volFracw_(volFracw)
{
}

Foam::waterEvaporationModel::waterEvaporationModel
(
    const dictionary& dict,
    const volScalarField& Sw,
    const volScalarField& eps,
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& pc,
    const volScalarField& volFracw
)
    :
    evapModel_(dict.lookup("evapModel")),
    kEvap_(dict.lookup("kEvap")),
    rho_(dict.lookupOrDefault("rho",dimensionedScalar("rho",dimDensity,1000))),
    Vm_(dict.lookupOrDefault("Vm",dimensionedScalar("Vm",dimVolume/dimMoles,1.8e-5))),
    as_(dict.lookupOrDefault("as",dimensionedScalar("as",dimArea/dimVolume,15))),
    lDiff_(dict.lookupOrDefault("lDiff",dimensionedScalar("lDiff",dimLength,2e-4))),
    Sw_(Sw),
    eps_(eps),
    T_(T),
    p_(p),
    pc_(pc),
    volFracw_(volFracw)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
Foam::waterEvaporationModel::~waterEvaporationModel()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * *//
Foam::tmp<Foam::volScalarField> Foam::waterEvaporationModel::pSatw()
{
    tmp<volScalarField> pSatw
    (
        new volScalarField
        (
            IOobject
            (
                "pSatw",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimPressure
        )
    );

    pSatw = dimensionedScalar(dimPressure, 611) * exp(17.27 * (T_ - dimensionedScalar(dimTemperature,298.15)) / (T_ - dimensionedScalar(dimTemperature,35.85))); // Clausius-Clapeyron equation

    return pSatw;
}

Foam::tmp<Foam::volScalarField> Foam::waterEvaporationModel::pVapw()
{
    tmp<volScalarField> pVapw
    (
        new volScalarField
        (
            IOobject
            (
                "pVapw",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimPressure
        )
    );

    pVapw.ref() = pSatw() * exp(-pc_*Vm_/(Foam::constant::physicoChemical::R * T_));

    return pVapw;
}

Foam::tmp<Foam::volScalarField> Foam::waterEvaporationModel::HEvap()
{
    tmp<volScalarField> HEvap
    (
        new volScalarField
        (
            IOobject
            (
                "HEvap",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimEnergy/dimMoles
        )
    );

    HEvap.ref() = dimensionedScalar(dimEnergy/dimMoles, 2.783e3)
        - dimensionedScalar(dimEnergy/dimMoles/dimTemperature, 5.166e-2) * T_
        - dimensionedScalar(dimEnergy/dimMoles/dimTemperature/dimTemperature, 3.622e-3) * T_ * T_;

    return HEvap;
}

Foam::tmp<Foam::volScalarField> Foam::waterEvaporationModel::evaporationRate()
{
    tmp<volScalarField> evaporationRate
    (
        new volScalarField
        (
            IOobject
            (
                "evaporationRate",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionedScalar("evaporationRate",dimDensity/dimTime, 0.0)
        )
    );

    if(evapModel_ == "firstOrder")
    {
        evaporationRate.ref() = kEvap_ * eps_ * (1-Sw_) * rho_ * (pVapw() - p_ * volFracw_) / p_;
    }
    else if(evapModel_ == "diffusive")
    {
        evaporationRate.ref() = binaryDiffusionCoeff::calculate("H2O","Air",p_,T_) * as_ * (pVapw() - p_ * volFracw_) 
            * dimensionedScalar(dimMass/dimMoles,1.8e-4) / constant::physicoChemical::R / T_
            / lDiff_;
    }
    else if(evapModel_ == "localEquilibrium")
    {
        const objectRegistry& or_ = p_.db();

        volScalarField& cw_ = or_.lookupObjectRef<volScalarField>("molConcentrationH2O");
        cw_ = pVapw() / (constant::physicoChemical::R * T_);
        const surfaceScalarField& phiConvg_ = or_.lookupObject<surfaceScalarField>("phiConvg");
        const surfaceScalarField& phiDiff_ = or_.lookupObject<surfaceScalarField>("molPhiH2O");
        evaporationRate.ref() = (eps_ * (1-Sw_) * fvc::ddt(cw_) + fvc::div(phiConvg_,cw_) + fvc::div(phiDiff_)) * dimensionedScalar(dimMass/dimMoles,1.8e-4);
    }
    else
    {
        FatalErrorIn("waterEvaporationModel::evaporationRate()")
            << "Unknown evaporation model " << evapModel_ << nl
            << "Valid options are: firstOrder, diffusive, localEquilibrium" << nl
            << abort(FatalError);
    }
    
    return evaporationRate;
}

Foam::tmp<Foam::volScalarField> Foam::waterEvaporationModel::QEvap()
{
    tmp<volScalarField> QEvap
    (
        new volScalarField
        (
            IOobject
            (
                "QEvap",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimEnergy/dimVolume/dimTime
        )
    );

    QEvap.ref() = evaporationRate() / dimensionedScalar(dimMass/dimMoles,1.8e4) // kg/mol
        * HEvap();

    return QEvap;
}