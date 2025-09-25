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
    Foam::binaryDiffusionCoeffField

Description
    A class for computing binary diffusion coefficients field for two species at a 
    certain point.

SourceFiles
    binaryDiffusionCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusionCoeffField_H
#define binaryDiffusionCoeffField_H

#include "fvMesh.H"
#include "volFields.H"
#include <vector>
#include <algorithm>
#include "gasCompositionSpace.H"
#include "dictionary.H"
#include "wordList.H"
#include "label.H"
#include "binaryDiffusionCoeff.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class binaryDiffusionCoeffField Declaration
\*---------------------------------------------------------------------------*/ 

class binaryDiffusionCoeffField
{
    protected:
        //- protected data

        //- reference to mesh
        const fvMesh& mesh_;

        //- reference to total pressure
        const volScalarField& p_;

        //- reference to temperature
        const volScalarField& T_;

        //- number of species
        label nSpecie_;

        //- list species names
        wordList specieNames_;

        //- array of field of binary diffusion coefficients
        std::vector<std::vector<volScalarField>> DAB_;
    
    public:
        //- constructors

        //- construct from nSpecie and specieNames
        binaryDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const label nSpecie,
            const wordList& specieNames
        );

        //- construct from gasCompositionSpace
        binaryDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const gasCompositionSpace& gasCompositionSpace
        );

        //- construct from dictionary
        binaryDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const dictionary& gasSpecieDict
        );

        //- destructor
        ~binaryDiffusionCoeffField();

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
        const volScalarField& DAB
        (
            word specieA,
            word specieB
        ) const
        {
            return DAB_[findSpecie(specieA)][findSpecie(specieB)];
        }

        //- return binary diffusion coefficient fields
        const std::vector<std::vector<volScalarField>>& DAB() const
        {
            return DAB_;
        }

        //-  find method
        label findSpecie(word specieName) const;



};

} // end namespace Foam

#endif