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
    Foam::binaryDiffusionCoeff

Description
    A class for computing binary diffusion coefficients for two species at a 
    certain point.

SourceFiles
    binaryDiffusionCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusionCoeff_H
#define binaryDiffusionCoeff_H

#include "scalar.H"
#include "word.H"
#include "OmegaTable.H"
#include "moleSpecificDataTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class binaryDiffusionCoeff Declaration
\*---------------------------------------------------------------------------*/

namespace binaryDiffusionCoeff
{
    extern volScalarField calculate
    (
        const word specieA,
        const word specieB,
        const volScalarField& p,
        const volScalarField& T
    );

}

} // end namespace Foam

#endif