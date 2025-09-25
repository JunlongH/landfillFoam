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

Class
    Foam::knudsenDiffusionCoeffField

Description
    A class for computing knudsen diffusion coefficients field in porous 
    media

SourceFiles
    knudsenDiffusionCoeffField.C
\*---------------------------------------------------------------------------*/

#ifndef knudsenDiffusionCoeffField_H
#define knudsenDiffusionCoeffField_H

#include <vector>
#include <algorithm>
#include "scalar.H"
#include "word.H"
#include "volFields.H"
#include "gasCompositionSpace.H"
#include "moleSpecificDataTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class knudsenDiffusionCoeffField Declaration
\*---------------------------------------------------------------------------*/ 

class knudsenDiffusionCoeffField
{
    private:
        //- private data

        //- reference to mesh
        const fvMesh& mesh_;

        //- reference to temperature field
        const volScalarField& T_;

        //- reference to gas permeability field
        const volScalarField& kg_;

        //- reference to effective permeability field
        const volScalarField& ke_;

        //- reference temperature
        const dimensionedScalar Tref_;

        //- reference molecular weight
        const scalar Mref_;

        //- reference air viscosity at 298.15 K
        const dimensionedScalar muRef_; // kg/(m.s)

        //- number of species
        label nSpecie_;

        //- list species names
        wordList specieNames_;

        //- array of field of knudsen diffusion coefficients
        std::vector<volScalarField> Dik_;

        //- type of klinkenberg coefficient model
        word klinkenbergType_;

    public:
        
        //- constructors

        //- construct from nSpecie and specieNames
        knudsenDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& T,
            const volScalarField& kg,
            const volScalarField& ke,
            const label nSpecie,
            const wordList& specieNames,
            const word klinkenbergType
        );

        //- construct from gasCompositionSpace
        knudsenDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& T,
            const volScalarField& kg,
            const volScalarField& ke,
            const gasCompositionSpace& gasCompSpace,
            const word klinkenbergType
        );

        //- construct from dictionary
        knudsenDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& T,
            const volScalarField& kg,
            const volScalarField& ke,
            const dictionary& gasSpecieDict
        );

        //- destructor
        ~knudsenDiffusionCoeffField();

        //- member functions

        //- return number of species
        label nSpecie() const
        {
            return nSpecie_;
        }

        //- return list of species names
        const wordList& specieNames() const
        {
            return specieNames_;
        }

        //- update knudsen diffusion coefficients
        void update();

        //- return knudsen diffusion coefficients
        const volScalarField& Dik(word specieName) const
        {
            return Dik_[findSpecie(specieName)];
        }

        //- return knudsen diffusion coefficients fields
        const std::vector<volScalarField>& Dik() const
        {
            return Dik_;
        }
        //- find specie index
        label findSpecie(word specieName) const;

        //- return klinkenberg coefficient
        tmp<volScalarField>  klinkenberg() const;

};


} // End namespace Foam


#endif