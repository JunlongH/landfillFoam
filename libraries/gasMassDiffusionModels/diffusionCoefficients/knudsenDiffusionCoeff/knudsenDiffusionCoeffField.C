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

#include "knudsenDiffusionCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//- construct from nSpecie and specieNames
Foam::knudsenDiffusionCoeffField::knudsenDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& T,
    const volScalarField& kg,
    const volScalarField& ke,
    const label nSpecie,
    const wordList& specieNames,
    const word klinkenbergType
)
    :
    mesh_(mesh),
    T_(T),
    kg_(kg),
    ke_(ke),
    Tref_(dimTemperature,298.15),
    Mref_(28.97),
    muRef_(dimDynamicViscosity,1.849e-5),
    nSpecie_(nSpecie),
    specieNames_(specieNames),
    klinkenbergType_(klinkenbergType)
{
    forAll(specieNames_,i)
    {
        Dik_.push_back
        (
            Foam::volScalarField
            (
                IOobject
                (
                    "Dik_" + specieNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionSet(0,2,-1,0,0,0,0)
            )
        );
    }
}

//- construct from gasCompositionSpace
Foam::knudsenDiffusionCoeffField::knudsenDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& T,
    const volScalarField& kg,
    const volScalarField& ke,
    const gasCompositionSpace& gasCompSpace,
    const word klinkenbergType
)
    :
    mesh_(mesh),
    T_(T),
    kg_(kg),
    ke_(ke),
    Tref_(dimTemperature,298.15),
    Mref_(28.97),
    muRef_(dimDynamicViscosity,1.849e-5),
    nSpecie_(gasCompSpace.nGasSpecie()),
    specieNames_(gasCompSpace.gasSpecieNames()),
    klinkenbergType_(klinkenbergType)
{
    forAll(specieNames_,i)
    {
        Dik_.push_back
        (
            Foam::volScalarField
            (
                IOobject
                (
                    "Dik_" + specieNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionSet(0,2,-1,0,0,0,0)
            )
        );
    }
}

//- construct from dictionary
Foam::knudsenDiffusionCoeffField::knudsenDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& T,
    const volScalarField& kg,
    const volScalarField& ke,
    const dictionary& gasSpecieDict
)
    :
    mesh_(mesh),
    T_(T),
    kg_(kg),
    ke_(ke),
    Tref_(dimTemperature,298.15),
    Mref_(28.97),
    muRef_(dimDynamicViscosity,1.849e-5),
    specieNames_(gasSpecieDict.lookup("gasSpecieNames")),
    klinkenbergType_(gasSpecieDict.lookupOrDefault<word>("klinkenbergType","Heid"))
{
    nSpecie_ = specieNames_.size();
    forAll(specieNames_,i)
    {
        Dik_.push_back
        (
            Foam::volScalarField
            (
                IOobject
                (
                    "Dik_" + specieNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionSet(0,2,-1,0,0,0,0)
            )
        );
    }
}

//- destructor
Foam::knudsenDiffusionCoeffField::~knudsenDiffusionCoeffField()
{}

// * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

void Foam::knudsenDiffusionCoeffField::update()
{
    forAll(specieNames_,i)
    {
        scalar molWeighti_ = moleSpecificDataTable::moleSpecificDataTable[specieNames_[i]][0];
        Dik_[i] = kg_ / muRef_ * pow(Mref_, 0.5) / pow(Tref_,0.5) *pow(T_,0.5)
            / pow(molWeighti_ , 0.5)
            * klinkenberg();
    }
}

Foam::tmp<Foam::volScalarField> Foam::knudsenDiffusionCoeffField::klinkenberg() const
{
    tmp<volScalarField> klinkenberg
    (   
        new volScalarField
        (
            IOobject
            (
                "klinkenbergCoefficient",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionSet(1,-1,-2,0,0,0,0)
        )
    );

    volScalarField::Boundary& klinkbf =  klinkenberg.ref().boundaryFieldRef();
    const volScalarField::Boundary& kebf =  ke_.boundaryField();
    const volScalarField::Boundary& kgbf =  kg_.boundaryField();

    if(klinkenbergType_=="Heid")
    {
        klinkenberg.ref().field()=0.11*pow(ke_,-0.39); // kl_ in m^2
        klinkbf = 0.11*pow(kebf,-0.39); 
    }
    else if(klinkenbergType_=="Reinecke")
    {
        klinkenberg.ref()=5.57*pow(kg_,-0.24); // kg_ in m^2
        klinkbf = 5.57*pow(kgbf,-0.24);
    }
    else if(klinkenbergType_=="Tanikawa")
    {
        klinkenberg.ref()=0.152*pow(ke_,-0.37); // kl_ in m^2
        klinkbf = 0.152*pow(kebf,-0.37);
    }
    else if(klinkenbergType_=="Jones")
    {
        klinkenberg.ref()=pow(ke_,-0.33); // kl_ in m^2
        klinkbf = pow(kebf,-0.33);
    }
    else
    {
        FatalErrorIn("knudsenDiffusionCoeffField::klinkenberg()")
            << "klinkenbergType_ should be either Heid, Reinecke, Tanikawa or Jones"
            << abort(FatalError);
    }
    return klinkenberg;
}

Foam::label Foam::knudsenDiffusionCoeffField::findSpecie(word specieName) const
{
    auto it = std::find(specieNames_.begin(),specieNames_.end(),specieName);
    if (it != specieNames_.end()) 
        {
        return std::distance(specieNames_.begin(), it);
        } 
    else 
        {
        // Element not found in the list
        FatalErrorIn("knudsenDiffusionCoeffField::findSpecie(const word& specieName) const")
            << "Specie " << specieName << " not found in the list of species names" << endl
            << abort(FatalError);
        return -1;
        }
}