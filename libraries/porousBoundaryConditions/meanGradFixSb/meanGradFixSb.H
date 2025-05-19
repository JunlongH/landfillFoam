/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::meanGradFixSb

Group
    porousBoundaryConditions

Description
    This boundary condition provides a gradient boundary condition for Sb of the
    two-phase flow solver. User provide a fixed value of Sb, and the boundary
    condition will calculate an gradient of the boundary and internal values.
          i.e. \ gradient_Sb = (Sb_b + Sb_i)/(2*distance_center_to_face)
    Different from normal fixed value boundary condition, this boundary
    condition for saturation uses a gradient in case of a dry or wet
    boundary causing a calculation error for capillary pressure and dpcdS.

Usage

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            meanGradFixSb;
        saturation      0.99;
    }
    \endverbatim

See also
    Foam::fixedGradient

SourceFiles
    meanGradFixSb.C

\*---------------------------------------------------------------------------*/

#ifndef meanGradFixSb_H
#define meanGradFixSb_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
        Class meanGradFixSbFvPatch Declaration
        \*---------------------------------------------------------------------------*/

    class meanGradFixSb
        :
        public fixedGradientFvPatchScalarField
    {
        // Private data

        // Protected data
    protected:
        //- Fixed value of saturation of the boundary condition
        scalar saturation_;


    public:

        //- Runtime type information
        TypeName("meanGradFixSb");


        // Constructors

        //- Construct from patch and internal field
        meanGradFixSb
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        meanGradFixSb
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //-  Construct by mapping given
        //  darcyGradPressure onto a new patch
        meanGradFixSb
        (
            const meanGradFixSb&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        meanGradFixSb
        (
            const meanGradFixSb&
        ) = delete;

        //- Construct as copy setting internal field reference
        meanGradFixSb
        (
            const meanGradFixSb&,
            const DimensionedField<scalar, volMesh>&
        ) ;
        
            //- Construct and return a clone setting internal field reference
            virtual tmp<fvPatchScalarField> clone
            (
                const DimensionedField<scalar, volMesh>& iF
            ) const
                {
                    return tmp<fvPatchScalarField>
                        (
                            new meanGradFixSb(*this, iF)
                        );
                }
        
        // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
