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
    Foam::landfillBiodegradation

Description
    Class of landfill biodegradation model.

SourceFiles
    landfillBiodegradation.C
    landfillBiodegradationNew.C
\*---------------------------------------------------------------------------*/

#ifndef landfillBiodegradation_H
#define landfillBiodegradation_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class landfillBiodegradation Declaration
\*---------------------------------------------------------------------------*/

class landfillBiodegradation
{
    //- protected data
    protected:

        word name_;
        dictionary biodegradProperties_;

        // reference to the substrate field
        volScalarField& substrate_;

        // Disallow copy constructor and assignment operator
        landfillBiodegradation
        (
            const landfillBiodegradation&
        );

        void operator=
        (
            const landfillBiodegradation&
        );

    public:

        //- Runtime type information
        TypeName("landfillBiodegradation");

        //- Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            landfillBiodegradation,
            dictionary,
            (
                const word& name,
                const dictionary& biodegradProperties,
                volScalarField& substrate
            ),
            (
                name,
                biodegradProperties,
                substrate
            )
        );

        //- Selectors
        static autoPtr<landfillBiodegradation> New
        (
            const word& name,
            const dictionary& biodegradProperties,
            volScalarField& substrate
        );

        //- Constructors
        landfillBiodegradation
        (
            const word& name,
            const dictionary& biodegradProperties,
            volScalarField& substrate
        );

        //- Destructor
        virtual ~landfillBiodegradation() {}

        //- Member Functions

        virtual volScalarField& biomassA() = 0;
        virtual volScalarField& biomassN() = 0;
        virtual volScalarField& inertMass() = 0;

        //- Return source term for substrate
        virtual tmp<volScalarField> sourceSubstrate() = 0;

        //virtual tmp<volScalarField> sourceSubstrateAerobic() = 0;

        //- Return source term for O2
        virtual tmp<volScalarField> sourceO2() = 0;

        //- Return source term for CO2
        virtual tmp<volScalarField> sourceCO2() = 0;

        //- Return source term for CH4
        virtual tmp<volScalarField> sourceCH4() = 0;

        //- Return source term for water
        virtual tmp<volScalarField> sourcew() = 0;

        //- Return heat source term caused by degradation
        virtual tmp<volScalarField> sourceHeat() = 0;

        //- Update member data
        virtual void update() = 0;



};
} // End namespace Foam




#endif