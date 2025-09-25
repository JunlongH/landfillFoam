#include "binaryDiffusionCoeff.H"

Foam::volScalarField Foam::binaryDiffusionCoeff::calculate
    (
        const word specieA,
        const word specieB,
        const volScalarField& p,
        const volScalarField& T
    )
    {
        scalar collisionDiameter_ = 0;
        scalar interactionT_ = 0;
        volScalarField omegaD_(T);
        volScalarField reducedT_(T);
        volScalarField DAB_
            (
                IOobject
                (
                    "DAB",
                    p.mesh().time().timeName(),
                    p.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                p.mesh(),
                dimensionedScalar
                    (
                        "DAB",
                        dimensionSet(0,2,-1,0,0,0,0),
                        0.0
                    )
            ); 
        collisionDiameter_= 0.5*(moleSpecificDataTable::moleSpecificDataTable[specieA][1] 
            + moleSpecificDataTable::moleSpecificDataTable[specieB][1]);

        interactionT_ = sqrt(moleSpecificDataTable::moleSpecificDataTable[specieA][2]
            * moleSpecificDataTable::moleSpecificDataTable[specieB][2]);

        reducedT_ = T/interactionT_;

        //- for every cell in the mesh interpolate
        forAll(reducedT_, cellI)
        {
            omegaD_[cellI] = OmegaTable::interpolate(OmegaTable::OmegaTable, reducedT_[cellI]);
        }

        forAll(reducedT_.boundaryField(), patchI)
        {
            forAll(reducedT_.boundaryField()[patchI], faceI)
            {
                omegaD_.boundaryFieldRef()[patchI][faceI] = OmegaTable::interpolate(OmegaTable::OmegaTable, reducedT_.boundaryField()[patchI][faceI]);
            }
        }

        DAB_.primitiveFieldRef()  = 0.0018583*pow(T.internalField(), 1.5)
            * sqrt(1/moleSpecificDataTable::moleSpecificDataTable[specieA][0] + 1/moleSpecificDataTable::moleSpecificDataTable[specieB][0])
            / (p.internalField()/101325*collisionDiameter_*collisionDiameter_*omegaD_.internalField()) // p at atm
            * 1e-4; // convert from cm2/s to m2/s
        
        DAB_.boundaryFieldRef() =  0.0018583*pow(T.boundaryField(), 1.5)
            * sqrt(1/moleSpecificDataTable::moleSpecificDataTable[specieA][0] + 1/moleSpecificDataTable::moleSpecificDataTable[specieB][0])
            / (p.boundaryField()/101325*collisionDiameter_*collisionDiameter_*omegaD_.boundaryField()) // p at atm
            * 1e-4; // convert from cm2/s to m2/s
        
        return DAB_;
    };