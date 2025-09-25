/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fluidGasSpecie.H"
#include "fixedValueFvPatchField.H"
#include "linear.H"
#include "fvc.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidGasSpecie::fluidGasSpecie
(
    const fvMesh& mesh,
    const dictionary& gasSpecieDict,
    const word& gasSpecieName,
    const volScalarField& p,
    const volScalarField& T
)
    :
    gasSpecie(mesh,gasSpecieDict,gasSpecieName),
    p_(p),
    T_(T),
    uSup_ 
    (
        IOobject
        (
            "uSup" + gasSpecieName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "uSup",
            dimensionSet(0,1,-1,0,0,0,0),
            vector::zero
        ),
        fvPatchField<vector>::calculatedType()
    ),
    molConcentration_
    (
        IOobject
        (
            "molConcentration" + gasSpecieName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    volFrac_
    (
        IOobject
        (
            "volFrac" + gasSpecieName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    massConcentration_
    (
        IOobject
        (
            "massConcentration" + gasSpecieName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        molConcentration_ * molWeight_ / 1000
    ),
    molPhiDiff_
    (
        IOobject
        (
            "molPhiDiff" + gasSpecieName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "molPhi",
            dimensionSet(0,0,-1,0,1,0,0),
            0.0
        )
    )
{

}

//- Fluid gas specie autoPtr
Foam::autoPtr<Foam::fluidGasSpecie> Foam::fluidGasSpecie::New
(
    const fvMesh& mesh,
    const dictionary& gasSpecieDict,
    const word& gasSpecieName,
    const volScalarField& p,
    const volScalarField& T
)
{
    return autoPtr<fluidGasSpecie>
    (
        new fluidGasSpecie(mesh, gasSpecieDict, gasSpecieName,p,T)
    );
}

//- Destructor
Foam::fluidGasSpecie::~fluidGasSpecie()
{}

// ************************************************************************* //

//- Member functions

//- Update velocity field
void Foam::fluidGasSpecie::updateMassConcentration()
{
    massConcentration_ = molConcentration_ * molWeight_ / 1000; // molweight in g/mol
}

//- Update volume fraction field
void Foam::fluidGasSpecie::updateVolFrac()
{
    volFrac_ = molConcentration_ * Foam::constant::physicoChemical::R * T_ / p_;
}


//- Update superficial velocity
void Foam::fluidGasSpecie::updateUSup()
{
    uSup_ = fvc::reconstruct(molPhiDiff_)/(molConcentration_+dimensionedScalar(dimensionSet(0,-3,0,0,1,0,0),VSMALL)); //$m\ s^{-1}$
}

void Foam::fluidGasSpecie::updateMolC()
{
    molConcentration_ = volFrac_ / Foam::constant::physicoChemical::R
        / T_ * p_; //$mol\ m^{-2}s^{-1}$
}

//- Update all fields
void Foam::fluidGasSpecie::updatebyMolC()
{
    updateMassConcentration();
    updateVolFrac();
    //updateUSup();
}