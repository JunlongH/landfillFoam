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
    Foam::gasCompositionSpce

Description
    A class for gas composition at any point of a field.

SourceFiles
    gasCompositionSpace.C
\*---------------------------------------------------------------------------*/

#ifndef gasCompositionSpace_H
#define gasCompositionSpace_H

#include "fluidGasSpecie.H"
#include "HashPtrTable.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    /*---------------------------------------------------------------------------*\
                   Class gasCompositionSpace Declaration
    \*---------------------------------------------------------------------------*/

    class gasCompositionSpace
    {
        //- protected data
        private:

            //- name of the gas composition space
            word name_;

            //- reference to mesh
            const fvMesh& mesh_;

            //- reference to total pressure
            volScalarField& p_;

            //- reference to temperature
            const volScalarField& T_;

            //- Is it necessary to store all the fields?
            //- total mole concentration
            volScalarField C_;

            //- total mass concentration
            volScalarField massC_;

            //- average molecular weight
            volScalarField molWeight_;

            //- average specific heat capacity
            volScalarField cp_;

            //- total flux
            surfaceScalarField phig_;

            //- total mass velocity
            //volVectorField U_;
            
            //- hash table of gas species
            HashPtrTable<fluidGasSpecie,word> gasSpecieTable_;

            //- number of gas species
            label nGasSpecie_;

            //- dictionary of gas species
            dictionary gasSpecieDict_;

            //- list of gas species names
            wordList gasSpecieNames_;



        //- public members
        public:

            //- Constructors

            gasCompositionSpace
            (
                const fvMesh& mesh,
                volScalarField& p,
                const volScalarField& T,
                const dictionary& gasSpecieDict
            );

            gasCompositionSpace
            (
                word name,
                const fvMesh& mesh,
                volScalarField& p,
                const volScalarField& T,
                const dictionary& gasSpecieDict
            );

            //- Destructor
            ~gasCompositionSpace();

            //- Member functions

            //- Return data member
            //- Return total mole concentration
            const volScalarField& molConcentration() const
            {
                return C_;
            }

            //- Return total mass concentration
            const volScalarField& massConcentration() const
            {
                return massC_;
            }

            //- Return average molecular weight
            const volScalarField& molWeight() const
            {
                return molWeight_;
            }

            //- Return average specific heat capacity
            const volScalarField& cp() const
            {
                return cp_;
            }

            //- Return total flux
            const surfaceScalarField& phig() const
            {
                return phig_;
            }

            //- Return velocity
            /*volVectorField& U()
            {
                return U_;
            }*/

            //- Return a reference of HashPtrTable
            HashPtrTable<fluidGasSpecie,word>& gasSpecieTable()
            {
                return gasSpecieTable_;
            }

            //- Return number of species
            label nGasSpecie() const
            {
                return nGasSpecie_;
            }

            //- Return list of species names
            const wordList& gasSpecieNames() const
            {
                return gasSpecieNames_;
            }

            //- Update concentration
            void updateMolConcentration();

            //- Update mass concentration
            void updateMassConcentration();

            //- Update pressure field
            void updatePressure();

            //- Update average molecular weight
            void updateMolWeight();

            //- Update average specific heat capacity
            void updateCp();
            
            //- Update velocity field
            void updateVelocity();

            //- Update all
            void update();

    
    };
} //End namespace Foam

//***************************************************************************//
#endif

//***************************************************************************//