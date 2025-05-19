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

\*---------------------------------------------------------------------------*/

#include "meanGradFixSb.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "linear.H"
#include "parRun.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "linear.H"
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "IOMRFZoneList.H"
#include "constants.H"

#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meanGradFixSb::meanGradFixSb
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    saturation_(0.0)
{
	fvPatchField<scalar>::operator=(patchInternalField());
	gradient() = 0.0;
}

Foam::meanGradFixSb::meanGradFixSb
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedGradientFvPatchScalarField(p, iF),
    saturation_(dict.lookupOrDefault<scalar>("saturation",0.0))
{}

Foam::meanGradFixSb::meanGradFixSb
(
	const meanGradFixSb& ptf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
	:
	fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    saturation_(ptf.saturation_)
{}


Foam::meanGradFixSb::meanGradFixSb
(
    const meanGradFixSb& ptf,
	const DimensionedField<scalar, volMesh>& iF
)
	:
	fixedGradientFvPatchScalarField(ptf, iF),
    saturation_(ptf.saturation_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meanGradFixSb::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& Cf = this->patch().Cf();
    const vectorField& Cn = this->patch().Cn();
    const vectorField& nf = this->patch().nf();
    const scalarField& Sb_i(this->patchInternalField());
    gradient() = (saturation_ - Sb_i)/2;
    gradient() /= mag((Cf - Cn) & nf);


    fixedGradientFvPatchScalarField::updateCoeffs();
}

void Foam::meanGradFixSb::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        meanGradFixSb
    );
}

// ************************************************************************* //
