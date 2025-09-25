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

#ifndef capillarityModel_H
#define capillarityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class capillarityModel Declaration
\*---------------------------------------------------------------------------*/

class capillarityModel
{

protected:

    // Protected data
    word name_;
    dictionary capillarityProperties_;

    // reference to the saturation field 
    const volScalarField& Sb_;

    // capillary pressure, derivative and capillary capacity
    volScalarField pc_,dpcdS_,Ch_;

    // Effective saturation field
    volScalarField Se_;

    // Private Member Functions

    //- Disallow copy construct
    capillarityModel(const capillarityModel&);

    //- Disallow default bitwise assignment
    void operator=(const capillarityModel&);

public:

    //- Runtime type information
    TypeName("capillarityModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        capillarityModel,
        dictionary,
        (
            const word& name,
            const dictionary& capillarityProperties,
            const volScalarField& Sb
        ),
        (name, capillarityProperties, Sb)
    );


    // Selectors

    //- Return a reference to the selected capillarity model
    static autoPtr<capillarityModel> New
    (
        const word& name,
        const dictionary& capillarityProperties,
        const volScalarField& Sb
    );


    // Constructors

    //- Construct from components
    capillarityModel
    (
        const word& name,
        const dictionary& capillarityProperties,
        const volScalarField& Sb
    );


    //- Destructor
    virtual ~capillarityModel() {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& capillarityProperties() const
        {
            return capillarityProperties_;
        }
      
    //- Return the capillary pressure as a function of Sb 
    tmp<volScalarField> pc() const
        {
            return pc_;
        }
    
    scalar pc(const label& celli)
        {
            return pc_[celli];
        }

    //- Return the derivative of capillary pressure as a function of Sb
    tmp<volScalarField> dpcdS() const
        {
            return dpcdS_;
        }
    scalar dpcdS(const label& celli)
        {
            return dpcdS_[celli];
        }
    //- Return effective saturation
    tmp<volScalarField> Se() const
        {
            return Se_;
        }

    //- Capillary capacity (Richards' model)
    tmp<volScalarField> Ch() const
        {
            return Ch_;
        };

    //- Capillary capacity (one cell)
    const scalar& Ch(const label& celli) const
        {
            return Ch_[celli];
        }

    //- To set specific value of Ch
    void setCh(const label& celli, const scalar& ChValue)
        {
            Ch_[celli] = ChValue;
        }
    
    //- Correct the capillary pressure
    virtual void correct() = 0; 
    virtual void correct(const label& celli) = 0;

    //- Compute capillary capacity from head pressure and return moisture (Richards' model)
    virtual volScalarField correctAndSb(volScalarField h) = 0;
    virtual scalar correctAndSb(const volScalarField& h, const label& celli) = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
