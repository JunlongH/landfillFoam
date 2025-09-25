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


Class
    Foam::atomicDiffusionVolume

Description
    A table of atomic diffusion volume for several elements

SourceFiles
    atomicDiffusionVolume.C

\*---------------------------------------------------------------------------*/

#ifndef atomicDiffusionVolume_H
#define atomicDiffusionVolume_H

#include "scalar.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class atomicDiffusionVolume Declaration
\*---------------------------------------------------------------------------*/

class atomicDiffusionVolumeTable
:
    public HashTable<scalar>
{

public:
  
    // Public types and data

        //- Structure to hold the element name and atomic diffusion volume pair
        struct atomicDiffusionVolume
        {
            char name[3];
            scalar diffusionVolume;
        };
        
        static const int nElements = 9; //- Pending to check

        //- Static table of the diffusion volume of several known elements
        static const atomicDiffusionVolume atomicDiffusionVolumes[nElements];


    // Constructors

        //- Construct from atomicDiffusionVolume_
        atomicDiffusionVolumeTable();
};


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// Atomic diffusion volume table for each element in the periodic table
extern atomicDiffusionVolumeTable atomicDiffusionVolumes;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
