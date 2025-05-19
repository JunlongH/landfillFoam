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
    Foam::methaneOxidation

Description
    Class for methane oxidation model, using Michaelis-Menten kinetics.

SourceFiles
    methaneOxidation.C
\*---------------------------------------------------------------------------*/

#ifndef methaneOxidation_H
#define methaneOxidation_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class methaneOxidation Declaration
\*---------------------------------------------------------------------------*/

class methaneOxidation
{
    //- protected data
    protected:
    
        dictionary methaneOxidationCoeffs_;

        const volScalarField& T_; // reference to temperature field
        const volScalarField& p_; // reference to pressure field
        const volScalarField& Sw_; // reference to water saturation field
        const volScalarField& molConcCH4_; // reference to methane mol concentration field
        const volScalarField& molConcO2_; // reference to oxygen mol concentration field


        scalar stoicCH4_; // stoichiometric coefficient of CH4
        scalar stoicO2_; // stoichiometric coefficient of oxygen degradating 1 unit of CH4
        scalar stoicCO2_; // stoichiometric coefficient of CO2 produced from 1 unit of CH4
        scalar stoicH2O_; // stoichiometric coefficient of H2O produced from 1 unit of CH4

        scalar Swmin_; // minimum water saturation
        scalar Swmax_; // maximum water saturation
        scalar Swfc_; // field capacity water saturation
        dimensionedScalar kmax_; // maximum specific rate of methane oxidation
        dimensionedScalar KmCH4_; // Michaelis-Menten constant for methane
        dimensionedScalar KmO2_; // Michaelis-Menten constant for oxygen
        dimensionedScalar enthalpyCH4_; // enthalpy of methane oxidation J/kg
        

        volScalarField kTemp_; // temperature correction factor for kmax
        volScalarField kSw_; // water saturation correction factor for kmax
        volScalarField kO2_; // oxygen concentration correction factor for kmax
        volScalarField kCH4_; // methane concentration correction factor for kmax

        void updatekTemp_(); // temperature correction factor for kmax
        void updatekSw_(); // water saturation correction factor for kmax
        void updatekO2_(); // oxygen concentration correction factor for kmax
        void updatekCH4_(); // methane concentration correction factor for kmax

    public:

        //- Constructor
        methaneOxidation
        (
            const dictionary& reactionProperties,
            const volScalarField& T,
            const volScalarField& p,
            const volScalarField& Sw
        );

        methaneOxidation
        (
            const dictionary& reactionProperties,
            const volScalarField& T,
            const volScalarField& p,
            const volScalarField& Sw,
            const volScalarField& molConcCH4,
            const volScalarField& molConcO2
        );
        //- Destructor
        ~methaneOxidation()
        {}

        //- Member Functions
        tmp<volScalarField> sourceCH4(); // mol/m3/s
        tmp<volScalarField> sourceO2(); // mol/m3/s
        tmp<volScalarField> sourceCO2(); // mol/m3/s
        tmp<volScalarField> sourcew(); // kg/m3/s
        tmp<volScalarField> sourceHeat(); // W/m3
        void update(); // update member data

};

} // End namespace Foam



#endif