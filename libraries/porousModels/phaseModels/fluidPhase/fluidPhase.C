/*---------------------------------------------------------------------------*\

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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

#include "fluidPhase.H"
#include "fixedValueFvPatchFields.H"
#include "linear.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhase::fluidPhase
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
    :  
    phase(mesh,transportProperties,phaseName),
    U_
    (
        IOobject
        (
            "U" + phaseName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
}


Foam::autoPtr<Foam::fluidPhase> Foam::fluidPhase::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
{
    return autoPtr<fluidPhase>
    (
        new fluidPhase(mesh, transportProperties, phaseName)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhase::~fluidPhase()
{}


// ************************************************************************* //
