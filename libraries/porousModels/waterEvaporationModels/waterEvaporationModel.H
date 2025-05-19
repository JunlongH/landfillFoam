/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::waterEvaporationModel

Description
    Calculation of water evaporation in porous media.

Source File
    waterEvaporationModel.C
\*---------------------------------------------------------------------------*/

#ifndef waterEvaporationModel_H
#define waterEvaporationModel_H

#include "volFields.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                             Class waterEvaporationModel Declaration
    \*---------------------------------------------------------------------------*/

class waterEvaporationModel
{
    private:
        //- Model type
        const word evapModel_;

        //- First order rate constant of water evaporation
        const dimensionedScalar kEvap_;

        //- Water density
        const dimensionedScalar rho_;

        //- Molar volume of water
        const dimensionedScalar Vm_;

        //- Specific area of water-gas interface
        const dimensionedScalar as_;

        //- Diffusion layer length of water
        const dimensionedScalar lDiff_;
        
        //- Reference to Water saturation field
        const volScalarField& Sw_;

        //- Reference to porosity field
        const volScalarField& eps_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to pressure field
        const volScalarField& p_;

        //- Reference to capillary pressure field
        const volScalarField& pc_;

        //- Reference to water vapor volumetric fraction field
        const volScalarField& volFracw_;
    
    public:
        //- Constructors
        //- Construct from components
        waterEvaporationModel
        (
            const word evapModel,
            const dimensionedScalar kEvap,
            const dimensionedScalar rho,
            const dimensionedScalar Vm,
            const dimensionedScalar as,
            const dimensionedScalar lDiff,
            const volScalarField& Sw,
            const volScalarField& eps,
            const volScalarField& T,
            const volScalarField& p,
            const volScalarField& pc,
            const volScalarField& volFracw
        );

        //- Construct from dictionary
        waterEvaporationModel
        (
            const dictionary& dict,
            const volScalarField& Sw,
            const volScalarField& eps,
            const volScalarField& T,
            const volScalarField& p,
            const volScalarField& pc,
            const volScalarField& volFracw
        );

        //- Destructor
        ~waterEvaporationModel();

        //- Member Functions
        //- Calculate water saturation pressure in free air
        tmp<volScalarField> pSatw();

        //- Calculate water saturation pressure in porous media
        tmp<volScalarField> pVapw();

        //- Calculate evaporation enthapy
        tmp<volScalarField> HEvap();

        //- Calculate water evaporation rate / kg/m^3/s
        tmp<volScalarField> evaporationRate();

        //- Calculate evaporation heat
        tmp<volScalarField> QEvap();
};
} // End namespace Foam

#endif