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

#include "binaryDiffusionCoeffField.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//- construct from nSpecie and specieNames
Foam::binaryDiffusionCoeffField::binaryDiffusionCoeffField
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
    //- resize the field
    // for each sub vector of DAB_ of the 2-D vector, resize it to nSpecie_
    for(int i = 0; i < nSpecie_; i++)
    {
        std::vector<Foam::volScalarField> row;
        for(int j = 0; j < nSpecie_; j++)
        {
            row.push_back
            (
                Foam::volScalarField
                (
                    IOobject
                    (
                        "DAB_" + specieNames_[i] + "_" + specieNames_[j],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                        (
                            "DAB_",
                            dimensionSet(0,2,-1,0,0,0,0),
                            0.0
                        ),
                    fvPatchField<scalar>::calculatedType()
                )
            );
        }
        DAB_.push_back(row);
    }
}

//- construct from gasCompositionSpace
Foam::binaryDiffusionCoeffField::binaryDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const gasCompositionSpace& gasComposition
)
    :
    mesh_(mesh),
    p_(p),
    T_(T),
    nSpecie_(gasComposition.nGasSpecie()),
    specieNames_(gasComposition.gasSpecieNames())
{
    for(int i = 0; i < nSpecie_; i++)
    {
        std::vector<Foam::volScalarField> row;
        for(int j = 0; j < nSpecie_; j++)
        {
            row.push_back
            (
                Foam::volScalarField
                (
                    IOobject
                    (
                        "DAB_" + specieNames_[i] + "_" + specieNames_[j],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                        (
                            "DAB_",
                            dimensionSet(0,2,-1,0,0,0,0),
                            0.0
                        ),
                    fvPatchField<scalar>::calculatedType()
                )
            );
        }
        DAB_.push_back(row);
    }
}

    //- construct from dictionary
Foam::binaryDiffusionCoeffField::binaryDiffusionCoeffField
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
    //- resize the field
    nSpecie_ = specieNames_.size();
    for(int i = 0; i < nSpecie_; i++)
    {
        std::vector<Foam::volScalarField> row;
        for(int j = 0; j < nSpecie_; j++)
        {
            row.push_back
            (
                Foam::volScalarField
                (
                    IOobject
                    (
                        "DAB_" + specieNames_[i] + "_" + specieNames_[j],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                        (
                            "DAB_",
                            dimensionSet(0,2,-1,0,0,0,0),
                            0.0
                        ),
                    fvPatchField<scalar>::calculatedType()
                )
            );
        }
        DAB_.push_back(row);
    }
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
Foam::binaryDiffusionCoeffField::~binaryDiffusionCoeffField()
{}


// * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //
Foam::label Foam::binaryDiffusionCoeffField::findSpecie(word specieName) const
{
    auto it = std::find(specieNames_.begin(),specieNames_.end(),specieName);
    if (it != specieNames_.end()) 
        {
        return std::distance(specieNames_.begin(), it);
        } 
    else 
        {
        // Element not found in the list
        FatalErrorIn("binaryDiffusionCoeffField::findSpecie(const word& specieName) const")
            << "Specie " << specieName << " not found in the list of species names" << endl
            << abort(FatalError);
        return -1;
        }
}

void Foam::binaryDiffusionCoeffField::update()
{
    //- update binary diffusion coefficients
    forAll(specieNames_,i)
    {
        forAll(specieNames_,j)
        {
            DAB_[i][j] = Foam::binaryDiffusionCoeff::calculate
            (
                specieNames_[i],
                specieNames_[j],
                p_,
                T_
            );
        }
    }
}
