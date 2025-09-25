/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/
#include "fickDiffusion.H"
#include "fvc.H"

Foam::List<Foam::tmp<Foam::surfaceScalarField>> Foam::gasMassDiffusion::fickDiffusion
(
    const fvMesh& mesh_,
    gasCompositionSpace& gasCompositionSpace_,
    const std::vector<volScalarField>& Di_
)
{
    const label N = gasCompositionSpace_.nGasSpecie();

    List<tmp<surfaceScalarField>> flux_(N);

    //- tmp for volFrac
    List<tmp<surfaceScalarField>> volFrac_(N); 
    //- tmp for molU at cell center
    List<tmp<surfaceScalarField>> molU_(N); //- molU at cell center
    //- tmp for gradient of molConcentration at cell center
    List<tmp<surfaceScalarField>> gradMolConcentration_(N);

    List<tmp<surfaceScalarField>> molPhi_(N);    


    // Calculate the flux
    forAll(gasCompositionSpace_.gasSpecieTable(), i)
    {
        word specieName_i = gasCompositionSpace_.gasSpecieNames()[i];

        molU_[i] = new surfaceScalarField
        (
            IOobject
            (
                "molU_" + specieName_i,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "molU_" + specieName_i,
                dimensionSet(0,-2,-1,0,1,0,0),
                0.0
            ),
            fvsPatchField<scalar>::calculatedType()
        );

        molPhi_[i] = new surfaceScalarField
        (
            IOobject
            (
                "molPhi_" + specieName_i,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "molPhi_" + specieName_i,
                dimensionSet(0,0,-1,0,1,0,0),
                0.0
            ),
            fvsPatchField<scalar>::calculatedType()
        );

        surfaceScalarField& molUi_ = molU_[i].ref();

        gradMolConcentration_[i] = fvc::snGrad(gasCompositionSpace_.gasSpecieTable()[specieName_i]->molConcentration());

        tmp<volScalarField> Ditmp_(Di_[i]);
        word specieName_ = gasCompositionSpace_.gasSpecieNames()[i];
        molUi_ = 
            - fvc::interpolate(Ditmp_,"diffusionCoeff") 
            * gradMolConcentration_[i];

        molPhi_[i].ref() = molU_[i].ref() * mesh_.magSf();
    }

    return molPhi_;   
}
