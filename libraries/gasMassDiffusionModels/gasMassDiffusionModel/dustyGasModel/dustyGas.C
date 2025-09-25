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

#include "dustyGas.H"
#include "fvc.H"

Foam::List<Foam::tmp<Foam::surfaceScalarField>> Foam::gasMassDiffusion::dustyGas
(
    const fvMesh& mesh_,
    gasCompositionSpace& gasCompositionSpace_,
    const std::vector<volScalarField>& DiK_,
    const std::vector<std::vector<volScalarField>>& Dij_
)
{
    const label N = gasCompositionSpace_.nGasSpecie();

    //- tmp for volFrac
    List<tmp<surfaceScalarField>> volFrac_(N); 
    //- tmp for molU at cell center
    List<tmp<surfaceScalarField>> molU_(N); //- molU at cell center
    //- tmp for gradient of molConcentration at cell center
    List<tmp<surfaceScalarField>> gradMolConcentration_(N);

    List<surfaceScalarField::Boundary*> molUbf;

    std::vector<surfaceScalarField> DiKf;

    std::vector<std::vector<surfaceScalarField>> Dijf;
    
    List<tmp<surfaceScalarField>> molPhi_(N);



    forAll(gasCompositionSpace_.gasSpecieTable(), i)
    {
        word specieName_i = gasCompositionSpace_.gasSpecieNames()[i];

        volFrac_[i] = fvc::interpolate(gasCompositionSpace_.gasSpecieTable()[specieName_i]->volFrac());

        gradMolConcentration_[i] = fvc::snGrad(gasCompositionSpace_.gasSpecieTable()[specieName_i]->molConcentration());

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

        //- Boundary field list
        molUbf.append(&(molU_[i].ref().boundaryFieldRef()));

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

        DiKf.push_back(fvc::interpolate(DiK_[i],"diffusionCoeff").ref());

        std::vector<surfaceScalarField> row;
        forAll(gasCompositionSpace_.gasSpecieTable(), j)
        {
            row.push_back
            (
                surfaceScalarField
                (
                    IOobject
                    (
                        "Dijf" + specieName_i + gasCompositionSpace_.gasSpecieNames()[j],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fvc::interpolate(Dij_[i][j],"diffusionCoeff")
                )
            );
        }
        Dijf.push_back(row);

    }   

    forAll(mesh_.magSf(),faceI)
    {
        //- Coefficient matrix A 
        scalarSquareMatrix A_(N,N);
        //- Constant vector b
        List<scalar> b_(N);
        //- Result vector molFlux
        List<scalar> U_(N);

        //- Construct the coefficient matrix A and constant vector b
        forAll(gasCompositionSpace_.gasSpecieTable(), i)
        {            
            forAll(gasCompositionSpace_.gasSpecieTable(), j)
            {
                if(i == j)
                {
                    forAll(gasCompositionSpace_.gasSpecieTable(), k)
                    {
                        if(k != i)
                        {
                            A_(i,j) += volFrac_[k]()[faceI] / Dijf[k][i][faceI];
                        }
                    }
                    A_(i,j) += 1/DiKf[i][faceI];
                }
                else
                {
                    A_(i,j) = - volFrac_[i]()[faceI] / Dijf[i][j][faceI];
                }
            }
            b_[i] = - gradMolConcentration_[i]()[faceI];
        }

        //- Solve the system of linear equations
        Foam::solve(U_, A_, b_);         

        //- Assign values to molU_
        forAll(gasCompositionSpace_.gasSpecieTable(), i)
        {
            molU_[i].ref()[faceI] = U_[i];
        }   
    }

    //- Solving equations for boundary surfaces
    forAll(gradMolConcentration_[0]().boundaryField(), patchI)
    {

        if(mesh_.boundaryMesh().types()[patchI] != "empty") //- How to treat empty boundary?
        {
            forAll(volFrac_[0]().boundaryField()[patchI], faceI)
            {


                //- Coefficient matrix A 
                scalarSquareMatrix A_(N,N);
                //- Constant vector b
                List<scalar> b_(N);
                //- Result vector molFlux
                List<scalar> U_(N);

                //- Construct the coefficient matrix A and constant vector b
                forAll(gasCompositionSpace_.gasSpecieTable(), i)
                {
                    forAll(gasCompositionSpace_.gasSpecieTable(), j)
                    {

                        if(i == j)
                        {
                            forAll(gasCompositionSpace_.gasSpecieTable(), k)
                            {
                                if(k != i)
                                {
                                    A_(i,j) += volFrac_[k].ref().boundaryField()[patchI][faceI] / Dijf[i][k].boundaryField()[patchI][faceI];
                                }
                            }
                            A_(i,j) += 1/DiKf[i].boundaryField()[patchI][faceI];
                        }
                        else
                        {
                            A_(i,j) = - volFrac_[i].ref().boundaryField()[patchI][faceI] / Dijf[i][j].boundaryField()[patchI][faceI];
                        }
                    }

                    b_[i] = - gradMolConcentration_[i]().boundaryField()[patchI][faceI];
                }

                

                //- Solve the system of linear equations
                Foam::solve(U_, A_, b_);



                

                //- Reconstruct molU_ field from U_
                forAll(gasCompositionSpace_.gasSpecieTable(), i)
                {
                    surfaceScalarField::Boundary& molUbf_i = *molUbf[i];
                    molUbf_i[patchI][faceI] = U_[i];
                }

            }
        }
    }

    forAll(gasCompositionSpace_.gasSpecieTable(), i)
    {
        molPhi_[i].ref() = molU_[i].ref() * mesh_.magSf();
    }

    return molPhi_;

}