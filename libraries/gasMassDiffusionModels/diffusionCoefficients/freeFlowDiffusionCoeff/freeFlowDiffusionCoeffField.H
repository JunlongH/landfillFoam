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
    Foam::freeFlowDiffusionCoeffField

Description
    A class for computing diffusion coefficients field within free flow. 

SourceFiles
    freeFlowDiffusionCoeffField.C
\*---------------------------------------------------------------------------*/

#ifndef freeFlowDiffusionCoeffField_H
#define freeFlowDiffusionCoeffField_H

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
             Class freeFlowDiffusionCoeffField Declaration
\*---------------------------------------------------------------------------*/ 

class freeFlowDiffusionCoeffField
{
    private:
        //- private data

        //- reference to mesh
        const fvMesh& mesh_;

        //- reference to total pressure field
        const volScalarField& p_;

        //- reference to temperature field
        const volScalarField& T_;

        //- number of species
        label nSpecie_;

        //- list species names
        wordList specieNames_;

        //- array of field of knudsen diffusion coefficients
        std::vector<volScalarField> Dia_;

    public:

        //- Constructors

        //- construct from nSpecie and specieNames
        freeFlowDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const label nSpecie,
            const wordList& specieNames
        );

        //- construct from gasCompositionSpace
        freeFlowDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const gasCompositionSpace& gasCompSpace
        );

        //- construct from dictionary
        freeFlowDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const dictionary& gasSpecieDict
        );

        //- destructor
        ~freeFlowDiffusionCoeffField();

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

        //- update binary diffusion coefficients
        void update();

        //- return binary diffusion coefficient of specieA and specieB
        const volScalarField& Dia
        (
            word specieName
        ) const
        {
            return Dia_[findSpecie(specieName)];
        }

        const std::vector<volScalarField>& Dia() const
        {
            return Dia_;
        }

        //-  find method
        label findSpecie(word specieName) const;

};
}



#endif