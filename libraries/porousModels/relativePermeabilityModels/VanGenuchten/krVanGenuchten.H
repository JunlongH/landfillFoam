/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

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

Class
    Foam::relativePermeabilityModels::krVanGenuchten

Description
    Standard Van Genuchten relative permeability model.

SourceFiles
    krVanGenuchten.C

\*---------------------------------------------------------------------------*/

#ifndef krVanGenuchten_H
#define krVanGenuchten_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krVanGenuchten Declaration
    \*---------------------------------------------------------------------------*/

class krVanGenuchten
:
    public relativePermeabilityModel
{
    //- Private data
    volScalarField Smin_;
    volScalarField Smax_;

    //- Van Genuchten coefficients
    dictionary krVanGenuchtenCoeffs_;
    volScalarField m_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

public:

    //- Runtime type information
    TypeName("VanGenuchten");

    // Constructors

    //- Construct from components
    krVanGenuchten
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );

    //- Destructor
    ~krVanGenuchten()
        {}

    // Member Functions

    //- Initiate parameters
    void init(volScalarField& Smin, volScalarField& Smax, volScalarField& m)
        {
            //- Assign the parameters
            Smin_ = Smin;
            Smax_ = Smax;
            m_ = m;
        }


    //- Correct the relative permeabilities
    void correct()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            kra_ = kramax_ * pow(1-Se_,0.5) * pow(1-pow(Se_,1/m_)+VSMALL,2*m_);
            krb_ = krbmax_ * pow(Se_,0.5) * pow(1-pow(1-pow(Se_,1/m_),m_),2);

            dkradS_ = - pow((1-pow(Se_,1/m_)+VSMALL),2*m_-1) * (-5*pow(Se_,1/m_+1)+4*pow(Se_,1/m_)+Se_);
            dkradS_ *= 1/(2*pow((1-Se_),0.5)*Se_);
            dkradS_ *=  1/(Smax_ - Smin_);
            dkradS_ *= kramax_;

            dkrbdS_ = 0.5 * (1-pow((1-pow(Se_,1/m_)),m_));
            dkrbdS_ *= ( 4 * pow(Se_,1/m_-0.5) * pow( (1-pow(Se_,1/m_)+VSMALL) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) / pow(Se_,0.5);
            dkrbdS_ *= 1/(Smax_ - Smin_);
            dkrbdS_ *= krbmax_;
        }
    void correcta()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            kra_ = kramax_ * pow(1-Se_,0.5) * pow(1-pow(Se_,1/m_),2*m_);
            dkradS_ = - pow((1-pow(Se_,1/m_)+VSMALL),2*m_-1) * (-5*pow(Se_,1/m_+1)+4*pow(Se_,1/m_)+Se_);
            dkradS_ *= 1/(2*pow((1-Se_),0.5)*Se_);
            dkradS_ *=  1/(Smax_ - Smin_);
            dkradS_ *= kramax_;
        }

    void correcta(const label& celli)
        {
            scalar Se = (Sb_[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            kra_[celli] = kramax_[celli] * pow(1-Se,0.5) * pow(1-pow(Se,1/m_[celli]),2*m_[celli]);
        }

    void correctb()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_ * pow(Se_,0.5) * pow(1-pow(1-pow(Se_,1/m_),m_),2);
            dkrbdS_ = 0.5 * (1-pow((1-pow(Se_,1/m_)),m_));
            dkrbdS_ *= ( 4 * pow(Se_,1/m_-0.5) * pow( (1-pow(Se_,1/m_)+VSMALL) , m_-1)) - ( pow((1-pow(Se_,1/m_)),m_) -1) / pow(Se_,0.5);
            dkrbdS_ *= 1/(Smax_ - Smin_);
            dkrbdS_ *= krbmax_;
        }
    void correctb(const label& celli)
        {
            scalar Se = (Sb_[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = krbmax_[celli] * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_[celli]),m_[celli]),2);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
