/*---------------------------------------------------------------------------*\

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
    Foam::incompressiblePhase

Description
    Incompressible phase model for the "porousMultiphaseFoam" toolbox.

SourceFiles
    incompressiblePhase.C

\*---------------------------------------------------------------------------*/

    #ifndef incompressiblePhase_H
    #define incompressiblePhase_H

    #include "fluidPhase.H"
    #include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                   Class incompressiblePhase Declaration
    \*---------------------------------------------------------------------------*/

    class incompressiblePhase
    :
    public fluidPhase
    {

        //- dynamic viscosity
        dimensionedScalar mu_;

        //- density
        dimensionedScalar rho_;

        //- phase fluxes
        surfaceScalarField phi_;

    public:

        // Constructors

        incompressiblePhase
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );

        // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<incompressiblePhase> New
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );

        //- Destructor
        virtual ~incompressiblePhase();

        // Member Functions
        const dimensionedScalar& mu() const
            {
                return mu_;
            }
        
        const dimensionedScalar& rho() const
            {
                return rho_;
            }
        
        const surfaceScalarField& phi() const
            {
                return phi_;
            }
        
        surfaceScalarField& phi()
            {
                return phi_;
            }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
