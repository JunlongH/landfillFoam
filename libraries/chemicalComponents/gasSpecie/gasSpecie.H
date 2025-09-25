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
    Foam::gasSpecie

Description
    Class of gas species, providing a class for the molecular weight and 
    molecular diffusion volume of the gas species.

SourceFiles
    gasSpecie.C
\*---------------------------------------------------------------------------*/

#ifndef gasSpecie_H
#define gasSpecie_H

#include "dictionary.H"
#include "fvMesh.H"
#include "specieElement.H"
#include "dimensionedType.H"

namespace Foam
{

    /*---------------------------------------------------------------------------* \
                             Class element Declaration
    \*---------------------------------------------------------------------------*/
    class gasSpecie
    {
        protected:

            //- reference to the mesh
            const fvMesh& mesh_;

            // Private data
            dictionary dict_;

            //- Name of element
            word name_;

            //- Molecular weight
            dimensionedScalar molWeight_;

            //- Heat capacity at constant pressure
            dimensionedScalar cp_;

        public:

            //- Constructors

            //- Construct from mesh, dictionary and gas specie name
            gasSpecie
            (
                const fvMesh& mesh,
                const dictionary& gasSpecieDict,
                const word& gasSpecieName
            );

            // - Constructor from mesh, specieComposition list and gas specie name
            gasSpecie
            (
                const fvMesh& mesh,
                const List<specieElement>& specieComposition,
                const word& gasSpecieName
            );

            //- Selectors
            static autoPtr<gasSpecie> New
            (
                const fvMesh& mesh,
                const dictionary& gasSpecieDict,
                const word& gasSpecieName
            );

            static autoPtr<gasSpecie> New
            (
                const fvMesh& mesh,
                const List<specieElement>& specieComposition,
                const word& gasSpecieName
            );

            // - Destructor
            virtual ~gasSpecie();

            // - Member functions

            //- Return name of the specie
            const word& name() const
            {
                return name_;
            }

            // - Return molecular weight of the specie
            const dimensionedScalar& molWeight() const
            {
                return molWeight_;
            }

            // - Return heat capacity at constant pressure of the specie
            const dimensionedScalar& cp() const
            {
                return cp_;
            }


            // - Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

    };

}


#endif