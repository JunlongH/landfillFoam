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
    Foam::monodDegradation

Description
    Subclass of landfill biodegradation model, dealing with monod degradation model.

SourceFiles
    monodDegradation.C
\*---------------------------------------------------------------------------*/
#ifndef monodDegradation_H
#define monodDegradation_H

#include "landfillBiodegradation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace landfillBiodegradationModels
{
/*---------------------------------------------------------------------------*\
                  Class monodDegradation Declaration
\*---------------------------------------------------------------------------*/

class monodDegradation
:
    public landfillBiodegradation
{
    //- Protected data
    protected:

        dictionary monodCoeffs_;

        dimensionedScalar biomassA0_; // initial aerobic biomass concentration
        dimensionedScalar biomassN0_; // initial anarobic biomass concentration
        volScalarField biomassA_; // aerobic biomass concentration
        volScalarField biomassN_; // anarobic biomass concentration
        volScalarField inertMass_; // inertMass after decay
        volScalarField& molConcO2_; // reference to oxygen concentration
        const volScalarField& T_; // reference to temperature field
        const volScalarField& p_; // reference to pressure field

        scalar stoicO2A_; // stoichiometric number of O2 for aerobic reaction
        scalar stoicO2N_; // stoichiometric number of O2 for anarobic reaction
        scalar stoicCO2A_; // stoichiometric number of CO2 for aerobic reaction
        scalar stoicCO2N_; // stoichiometric number of CO2 for anarobic reaction
        scalar stoicCH4A_; // stoichiometric number of CH4 for aerobic reaction
        scalar stoicCH4N_; // stoichiometric number of CH4 for anarobic reaction
        scalar stoicH2OA_; // stoichiometric number of H2O for aerobic reaction
        scalar stoicH2ON_; // stoichiometric number of H2O for anarobic reaction



        dimensionedScalar kMaxA_; // maximum aerobic growth rate
        dimensionedScalar kMaxN_; // maximum anarobic growth rate
        dimensionedScalar KSA_; // half saturation constant for substrate for aerobic growth
        dimensionedScalar KSN_; // half saturation constant for substrate for anarobic growth
        dimensionedScalar KO2_; // half saturation constant for oxygen
        dimensionedScalar YSA_; // yield coefficient for aerobic growth
        dimensionedScalar YSN_; // yield coefficient for anarobic growth
        dimensionedScalar productFactorA_; // Yield factor for aerobic gas production
        dimensionedScalar productFactorN_; // Yield factor for anarobic gas production
        dimensionedScalar molWeightSubstrate_; // molecular weight of substrate
        dimensionedScalar molConcO2Anaerobic_; // anaerobic oxygen concentration
        dimensionedScalar decayConstant_; // decay constant for biomass

        dimensionedScalar enthalpyA_; // enthalpy of aerobic degradation
        dimensionedScalar enthalpyN_; // enthalpy of anarobic degradation

        dimensionedScalar TminA_; // minimum temperature for aerobic growth
        dimensionedScalar TmaxA_; // maximum temperature for aerobic growth
        dimensionedScalar ToptA_; // optimum temperature for aerobic growth
        dimensionedScalar TminN_; // minimum temperature for anarobic growth
        dimensionedScalar TmaxN_; // maximum temperature for anarobic growth
        dimensionedScalar ToptN_; // optimum temperature for anarobic growth

        volScalarField kTempA_; // temperature correction factor for aerobic growth
        volScalarField kTempN_; // temperature correction factor for anarobic growth

        volScalarField growthA_; // growth rate for aerobic biomass
        volScalarField growthN_; // growth rate for anarobic biomass
        volScalarField decayA_; // decay rate for aerobic biomass
        volScalarField decayN_; // decay rate for anarobic biomass

        void updatekTempA_(); // temperature correction factor for aerobic growth 
        void updatekTempN_(); // temperature correction factor for anarobic growth
        tmp<volScalarField> kSat_(volScalarField& s, dimensionedScalar& k); // saturation correction factor

        void updategrowthA(); // growth rate for aerobic biomass
        void updategrowthN(); // growth rate for anarobic biomass
        void updatedecayA(); // decay rate for aerobic biomass
        void updatedecayN(); // decay rate for anarobic biomass


    
    public:

        //- Runtime type information
        TypeName("monod");

        //- Constructor
        monodDegradation
        (
            const word& name,
            const dictionary& biodegradProperties,
            volScalarField& substrate
        );

        //- Destructor
        ~monodDegradation()
        {}

        //- Member Functions
        //- return the biomass concentration
        volScalarField& biomassA() {return biomassA_;}
        volScalarField& biomassN() {return biomassN_;}
        volScalarField& inertMass() {return inertMass_;}

        tmp<volScalarField> sourceSubstrate(); // kg/m3/s

        //tmp<volScalarField> sourceSubstrateAerobic(); // kg/m3/s

        tmp<volScalarField> sourceO2(); // mol/m3/s

        tmp<volScalarField> sourceCO2(); // mol/m3/s

        tmp<volScalarField> sourceCH4(); // mol/m3/s

        tmp<volScalarField> sourcew(); // kg/m3/s

        tmp<volScalarField> sourceHeat(); // W/m3

        void update(); // Solving biomass and substrate concentration

};
}// End namespace landfillBiodegradationModels
}// End namespace Foam



#endif