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
    Foam::porousBinaryDiffusionCoeffField

Description
    A class for computing binary diffusion coefficients field in porous media
    for two species at a certain point.

SourceFiles
    porousBinaryDiffusionCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBinaryDiffusionCoeffField_H
#define porousBinaryDiffusionCoeffField_H

#include "binaryDiffusionCoeffField.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class porousBinaryDiffusionCoeffField Declaration
\*---------------------------------------------------------------------------*/ 

class porousBinaryDiffusionCoeffField
:
    public binaryDiffusionCoeffField
{
    protected:
        //- protected data

        //- reference to porosity field
        const volScalarField& eps_;

        //- reference to gas saturation field
        const volScalarField& Sg_;

        //- reference to porous media tortuosity power calculating effective diffusion coefficient
        const volScalarField m_;

    
    public:
        //- Constructors

        //- Construct from porosity field and gas saturation field
        porousBinaryDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const volScalarField& eps,
            const volScalarField& Sg,
            const label nSpecie,
            const wordList& specieNames
        );

        porousBinaryDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const volScalarField& eps,
            const volScalarField& Sg,
            const gasCompositionSpace& gasCompositionSpace
        );

        porousBinaryDiffusionCoeffField
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volScalarField& T,
            const volScalarField& eps,
            const volScalarField& Sg,
            const dictionary& gasSpecieDict
        );

        //- Destructor
        ~porousBinaryDiffusionCoeffField();

        //- Member Functions      
        //- Update
        void update();

};


} // End namespace Foam

#endif