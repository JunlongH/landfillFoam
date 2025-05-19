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

\*---------------------------------------------------------------------------*/

#ifndef pcVanGenuchten_H
#define pcVanGenuchten_H

#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
                    Class pcVanGenuchten Declaration
\*---------------------------------------------------------------------------*/

class pcVanGenuchten
:
    public capillarityModel
{
    // pcVan Genuchten coefficients
    dictionary pcVanGenuchtenCoeffs_;
    volScalarField Smin_,Smax_;
    volScalarField m_,n_,alpha_;
    volScalarField pc0_;
        
public:

    //- Runtime type information
    TypeName("VanGenuchten");

    // Constructors

    //- Construct from components
    pcVanGenuchten
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& Sb
    );

    //- Destructor
    ~pcVanGenuchten()
        {}

    // Member Functions

    //- Return reference of parameters
    volScalarField& Smin()
        {
            return Smin_;
        }
    
    volScalarField& Smax()
        {
            return Smax_;
        }

    volScalarField& m()
        {
            return m_;
        }
    
    volScalarField& n()
        {
            return n_;
        }
    
    volScalarField& alpha()
        {
            return alpha_;
        }
    
    volScalarField& pc0()
        {
            return pc0_;
        }

    //- Initiate parameters
    void init(volScalarField& Smin, volScalarField& Smax, volScalarField& m, volScalarField& n, volScalarField& alpha, volScalarField& pc0)
        {
            Smin_ = Smin;
            Smax_ = Smax;
            m_ = m;
            n_ = n;
            alpha_ = alpha;
            pc0_ = pc0;
        }

    //- Correct the capillary pressure
    void correct()
        {
            Se_ = (Sb_-Smin_)/(Smax_-Smin_);
            pc_ = pc0_ * pow(pow(Se_,-1/m_)-1,1/n_);
            dpcdS_ = - (1/(n_*m_)) * (pc0_/(Smax_-Smin_)) * pow(pow(Se_+VSMALL,-1/m_)-1,(1/n_)-1) * pow(Se_+VSMALL,-(1+m_)/m_);

        }
    
    void correct(const label& celli)
        {
            scalar Se = (Sb_[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            
            pc_[celli] = pc0_[celli] * pow(pow(Se,-1/m_[celli])-1,1/n_[celli]);
            dpcdS_[celli] = - (1/(n_[celli]*m_[celli])) * (pc0_[celli]/(Smax_[celli]-Smin_[celli])) * pow(pow(Se+VSMALL,-1/m_[celli])-1,(1/n_[celli])-1) * pow(Se+VSMALL,-(1+m_[celli])/m_[celli]);
        }

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    volScalarField correctAndSb(volScalarField h)
        {
            //- Compute Saturation from! head pressure
            h.dimensions().reset(dimless);
            Se_ =  neg(h+VSMALL) *pow(1.0+pow(alpha_*mag(h),n_),-m_) + pos(h+VSMALL)*1;
            h.dimensions().reset(dimLength);

            //- Update capillary capacity
            Ch_.dimensions().reset(dimless);
            Ch_ = alpha_*m_*(Smax_-Smin_)/(1.0-m_)*pow(Se_,1.0/m_)*pow(1.0-pow(Se_,1.0/m_),m_) ;
            Ch_.dimensions().reset(dimless/dimLength);

            //- Update saturation
            volScalarField Stmp_(Smin_+(Smax_-Smin_)*Se_);
            
            return Stmp_;
        }

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(const volScalarField& h, const label& celli)
        {
            //- Compute Saturation from! head pressure
            scalar Se =  neg(h[celli]) *pow(1.0+pow(alpha_[celli]*mag(h[celli]),n_[celli]),-m_[celli]) + pos(h[celli]+VSMALL)*(1);

            //- Update capillary capacity
            Ch_[celli] = alpha_[celli]*m_[celli]*(Smax_[celli]-Smin_[celli])/(1.0-m_[celli])*pow(Se,1.0/m_[celli])*pow(1.0-pow(Se,1.0/m_[celli]),m_[celli]) ;

            //- Return moisture value
            return Smin_[celli]+(Smax_[celli]-Smin_[celli])*Se;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
