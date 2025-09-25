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

#include "incompressiblePhase.H"
#include "fixedValueFvPatchFields.H"
#include "linear.H"
#include "GeometricField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressiblePhase::incompressiblePhase
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
:
    fluidPhase(mesh,transportProperties,phaseName),
    mu_(dict_.lookup("mu")),
    rho_(dict_.lookup("rho")),
    phi_
    (
        IOobject
        (
            "phi" + phaseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,3,-1,0,0,0,0), 0.0)
    )
{   
}


Foam::autoPtr<Foam::incompressiblePhase> Foam::incompressiblePhase::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
{
    return autoPtr<incompressiblePhase>
    (
        new incompressiblePhase(mesh, transportProperties, phaseName)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressiblePhase::~incompressiblePhase()
{}


// ************************************************************************* //
