/*---------------------------------------------------------------------------*\
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

#include "element.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::element::element
(
    const fvMesh& mesh,
    const dictionary& elements,
    const word& elementName
)
    :
    mesh_(mesh),
    dict_(),
    name_(elementName)
{
    dict_ = elements.subDict("element." + elementName);
    atomicWeight_ = dict_.lookupOrDefault<scalar>("atomicWeight", 0);
    atomicDiffusionVolume_ = dict_.lookupOrDefault<scalar>("atomicDiffusionVolume", 0);
}
// - Element autoPtr
Foam::autoPtr<Foam::element> Foam::element::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& elementName
)
{
    return autoPtr<element>
        (
            new element(mesh, transportProperties, "element." + elementName)
        );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::element::~element()
{}

