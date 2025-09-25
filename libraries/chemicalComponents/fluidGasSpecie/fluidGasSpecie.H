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
    Foam::fluidGasSpecie

Description
    Class of gas species fluid, providing a class for fields of the gas 
    species.

SourceFiles
    fluidGasSpecie.C
\*---------------------------------------------------------------------------*/

#ifndef fluidGasSpecie_H
#define fluidGasSpecie_H

#include "gasSpecie.H"
#include "dimensionedType.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                   Class fluidGasSpecie Declaration
    \*---------------------------------------------------------------------------*/

    class fluidGasSpecie
    :
        public gasSpecie
    {
        //- protected data
        protected:

            //- Reference to pressure field
            const volScalarField& p_;

            //- Reference to temperature field
            const volScalarField& T_;

            //- specie superficial velocity
            volVectorField uSup_;

            //- specie mole concentration
            volScalarField molConcentration_; // mol/m3

            //- specie volume fraction
            volScalarField volFrac_; 

            //- specie mass concentration
            volScalarField massConcentration_; // kg/m3

            //- specie face mole flux
            //- This model calculates basing on face mole flux
            surfaceScalarField molPhiDiff_; // mol/s 

        //- public members
        public:

            //- Constructors

            fluidGasSpecie
            (
                const fvMesh& mesh,
                const dictionary& gasSpecieDict,
                const word& specieName,
                const volScalarField& p,
                const volScalarField& T
            );

            //- Selectors

            //- Return a reference to the specie
            static autoPtr<fluidGasSpecie> New
            (
                const fvMesh& mesh,
                const dictionary& gasSpecieDict,
                const word& specieName,
                const volScalarField& p,
                const volScalarField& T
            );

            //- Destructor
            virtual ~fluidGasSpecie();

            //- Member Functions

            //- Return velocity field
            const volVectorField& uSup() const
                {
                    return uSup_;
                }
            
            volVectorField& uSup()
                {
                    return uSup_;
                }
            
            //- Return molecular concentration field
            const volScalarField& molConcentration() const
                {
                    return molConcentration_;
                }
            
            volScalarField& molConcentration()
                {
                    return molConcentration_;
                }
            
            //- Return volume fraction field
            const volScalarField& volFrac() const
                {
                    return volFrac_;
                }
            
            volScalarField& volFrac()
                {
                    return volFrac_;
                }
            
            //- Return mass concentration field
            const volScalarField& massConcentration() const
                {
                    return massConcentration_;
                }

            volScalarField& massConcentration()
                {
                    return massConcentration_;
                }

            //- Return face mole flux field
            const surfaceScalarField& molPhiDiff() const
                {
                    return molPhiDiff_;
                }
            
            surfaceScalarField& molPhiDiff()
                {
                    return molPhiDiff_;
                }
            

            
            //- Update mass concentration field
            void updateMassConcentration();

            //- Update volume fraction field
            void updateVolFrac();

            //- Update superficial velocity field
            void updateUSup();

            //- Update mole concentration field
            void updateMolC();

            //- Update all fields
            void updatebyMolC();

    };

} // End namespace Foam

// ***************************************************************************//

#endif

//***************************************************************************//