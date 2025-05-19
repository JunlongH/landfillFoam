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
    Foam::autoExplicitRelaxFactor

Description
    Calculating factor basing on angle between two iterations.

\*---------------------------------------------------------------------------*/

#ifndef autoExplicitRelaxFactor_H
#define autoExplicitRelaxFactor_H

#include "fvMesh.H"
#include "volFieldsFwd.H"

namespace Foam
{
    template <typename T>
    scalar autoExplicitRelaxFactor(const VolField<T>& f, const VolField<T>& deltaf) 
    {
        scalar fRelaxFactorTimes = 1.0;

        Field<T> tdeltaF = deltaf.internalField();
        Field<T> tdeltaFPrevIter = deltaf.prevIter().internalField();
        Field<T> deltafProduct = tdeltaF * tdeltaFPrevIter;

        Field<T> deltafAve = pow(tdeltaF,2);
        Field<T> deltafPrevIterAve = pow(tdeltaFPrevIter,2);

        scalar cosDelta = (gSum(deltafProduct))
            /(sqrt(gSum(deltafAve))*sqrt(gSum(deltafPrevIterAve))+VSMALL);

        if(cosDelta > 0.7071)
        {
            fRelaxFactorTimes = 1.4142;
        }
        else if(cosDelta <= 0.7071 && cosDelta > 0.0)
        {
            fRelaxFactorTimes = 1.0;
        }
        else
        {
            fRelaxFactorTimes = 0.7071;
        }

        return fRelaxFactorTimes;
    }

    template<>
    scalar autoExplicitRelaxFactor(const VolField<vector>& f, const VolField<vector>& deltaf)
    {
        scalar fRelaxFactorTimes = 1.0;

        Field<scalar> tdeltaF = mag(deltaf.internalField());
        Field<scalar> tdeltaFPrevIter = mag(deltaf.prevIter().internalField());
        Field<scalar> deltafProduct = tdeltaF * tdeltaFPrevIter;

        Field<scalar> deltafAve = pow(tdeltaF,2);
        Field<scalar> deltafPrevIterAve = pow(tdeltaFPrevIter,2);

        scalar cosDelta = (gSum(deltafProduct))
            /(sqrt(gSum(deltafAve))*sqrt(gSum(deltafPrevIterAve))+VSMALL);

        if(cosDelta > 0.7071)
        {
            fRelaxFactorTimes = 1.4142;
        }
        else if(cosDelta <= 0.7071 && cosDelta > 0.0)
        {
            fRelaxFactorTimes = 1.0;
        }
        else
        {
            fRelaxFactorTimes = 0.7071;
        }

        return fRelaxFactorTimes;
    }


}

#endif