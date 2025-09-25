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

#include "atomicDiffusionVolumes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::atomicDiffusionVolumeTable::atomicDiffusionVolume
Foam::atomicDiffusionVolumeTable::atomicDiffusionVolumes[atomicDiffusionVolumeTable::nElements] =
{
    {"H",    2.31},
    {"He",   2.67},
    {"C",   15.9},
    {"N",   4.54},
    {"O",   6.11},
    {"S",   22.9},
    {"Cl",  21},
    {"Br",  21.9},
    {"I",   29.8},
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atomicDiffusionVolumeTable::atomicDiffusionVolumeTable()
{
    for (int i=0; i<nElements; i++)
    {
        insert(word(atomicDiffusionVolumes[i].name), atomicDiffusionVolumes[i].diffusionVolume);
    }
}


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

Foam::atomicDiffusionVolumeTable Foam::atomicDiffusionVolumes;


// ************************************************************************* //
