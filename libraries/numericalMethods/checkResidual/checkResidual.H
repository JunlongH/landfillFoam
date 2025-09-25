/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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


Function
    Foam::checkResidual

Description
    Checking residuals for convergence.
    
\*---------------------------------------------------------------------------*/

#ifndef checkResidual_H
#define checkResidual_H

#include "fvMesh.H"

namespace Foam
{
    template <typename T>
    bool checkResidual(const T& f, const T& deltaf, const scalar tolerance)
    {
        tmp<volScalarField> tmagDeltaF = mag(deltaf);
        tmp<volScalarField> tmagF = mag(f);

        scalar fMax = gMax(tmagF.ref());
        scalar deltafMax = gMax(tmagDeltaF.ref());

        scalar fResidual = deltafMax/(fMax+VSMALL);

        return (fResidual < tolerance); // return true if converged
    }
}

#endif