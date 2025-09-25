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

Class
    Foam::element

Description
    Class of chemical elements

SourceFiles
    element.C

\*---------------------------------------------------------------------------*/

#ifndef element_H
#define element_H

#include "dictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------* \
                             Class element Declaration
    \*---------------------------------------------------------------------------*/

    class element
    {

    protected:

        //- reference to the mesh
        const fvMesh& mesh_;

        // Private data
        dictionary dict_;

        //- Name of element
        word name_;

        //- Atomic weight
        scalar atomicWeight_;

        //- Atomic diffusion volume
        scalar atomicDiffusionVolume_;

    public:

        // Constructors

        element
        (
            const fvMesh& mesh,
            const dictionary& elements,
            const word& elementName
        );

        // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<element> New
        (
            const fvMesh& mesh,
            const dictionary& elements,
            const word& elementName
        );

        //- Destructor
        virtual ~element();

        // Member Functions
        const word& name() const
        {
            return name_;
        }

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        const scalar& atomicWeight() const
		{
			return atomicWeight_;
		}

        const scalar& atomicDiffusionVolume() const
        {
            return atomicDiffusionVolume_;
        }

    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
