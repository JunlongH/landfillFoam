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

#include "freeFlowDiffusionCoeffField.H"
#include "binaryDiffusionCoeff.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//- construct from nSpecie and specieNames
Foam::freeFlowDiffusionCoeffField::freeFlowDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const label nSpecie,
    const wordList& specieNames
)
    :
    mesh_(mesh),
    p_(p),
    T_(T),
    nSpecie_(nSpecie),
    specieNames_(specieNames)
{
    for(int i = 0; i < nSpecie_; i++)
        Dia_[i] = Foam::volScalarField
        (
            IOobject
            (
                "Dia_" + specieNames_[i],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionSet(0,2,-1,0,0,0,0)
        );
}

//- construct from gasCompositionSpace
Foam::freeFlowDiffusionCoeffField::freeFlowDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const gasCompositionSpace& gasCompSpace
)
    :
    mesh_(mesh),
    p_(p),
    T_(T),
    nSpecie_(gasCompSpace.nGasSpecie()),
    specieNames_(gasCompSpace.gasSpecieNames())
{
    for(int i = 0; i < nSpecie_; i++)
        Dia_[i] = Foam::volScalarField
        (
            IOobject
            (
                "Dia_" + specieNames_[i],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionSet(0,2,-1,0,0,0,0)
        );
}

//- construct from dictionary
Foam::freeFlowDiffusionCoeffField::freeFlowDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const dictionary& gasSpecieDict
)
    :
    mesh_(mesh),
    p_(p),
    T_(T),
    specieNames_(gasSpecieDict.lookup("specieNames"))
{
    nSpecie_ = specieNames_.size();
    for(int i = 0; i < nSpecie_; i++)
        Dia_[i] = Foam::volScalarField
        (
            IOobject
            (
                "Dia_" + specieNames_[i],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionSet(0,2,-1,0,0,0,0)
        );
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
Foam::freeFlowDiffusionCoeffField::~freeFlowDiffusionCoeffField()
{}

// * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //
Foam::label Foam::freeFlowDiffusionCoeffField::findSpecie(word specieName) const
{
    auto it = std::find(specieNames_.begin(),specieNames_.end(),specieName);
    if (it != specieNames_.end()) 
        {
        return std::distance(specieNames_.begin(), it);
        } 
    else 
        {
        // Element not found in the list
        FatalErrorIn("freeFlowDiffusionCoeffField::findSpecie(const word& specieName) const")
            << "Specie " << specieName << " not found in the list of species names" << endl
            << abort(FatalError);
        return -1;
        }
}

//- update binary diffusion coefficients
void Foam::freeFlowDiffusionCoeffField::update()
{
    forAll(specieNames_,i)
    {
        Dia_[i] = Foam::binaryDiffusionCoeff::calculate
        (
            specieNames_[i],
            "air",
            p_,
            T_
        );
    }
}
