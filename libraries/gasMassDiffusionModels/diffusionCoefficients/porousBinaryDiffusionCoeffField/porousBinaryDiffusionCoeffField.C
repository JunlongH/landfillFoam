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

#include "porousBinaryDiffusionCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- construct from nSpecie and specieNames
Foam::porousBinaryDiffusionCoeffField::porousBinaryDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const volScalarField& eps,
    const volScalarField& Sg,
    const label nSpecie,
    const wordList& specieNames
)
    :
    binaryDiffusionCoeffField(mesh, p, T, nSpecie, specieNames),
    eps_(eps),
    Sg_(Sg),
    m_(
        IOobject
        (
            "tortuosityParameter",
            p.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 4/3)
    )
{
}

Foam::porousBinaryDiffusionCoeffField::porousBinaryDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const volScalarField& eps,
    const volScalarField& Sg,
    const gasCompositionSpace& gasCompositionSpace
)
    :
    binaryDiffusionCoeffField(mesh, p, T, gasCompositionSpace),
    eps_(eps),
    Sg_(Sg),
    m_(
        IOobject
        (
            "tortuosityParameter",
            p.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 4/3)
    )
{
}

Foam::porousBinaryDiffusionCoeffField::porousBinaryDiffusionCoeffField
(
    const fvMesh& mesh,
    const volScalarField& p,
    const volScalarField& T,
    const volScalarField& eps,
    const volScalarField& Sg,
    const dictionary& gasSpecieDict
)
    :
    binaryDiffusionCoeffField(mesh, p, T, gasSpecieDict),
    eps_(eps),
    Sg_(Sg),
    m_(
        IOobject
        (
            "tortuosityParameter",
            p.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 4/3)
    )
{
}

//- Destructor
Foam::porousBinaryDiffusionCoeffField::~porousBinaryDiffusionCoeffField()
{
}

//- Update
void Foam::porousBinaryDiffusionCoeffField::update()
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
            DAB_[i][j] = DAB_[i][j]* pow(eps_,m_)
                * pow(Sg_,10/3); //Millington and Quirk (1961)
        }
    }
}




