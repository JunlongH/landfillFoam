//- Thermal diffusivity
DT = DTs*DTw*DTg/(eps*Sg*DTs*DTw+(1-eps)*DTw*DTg+eps*Sw*DTg*DTs);

//- Initiate species relatied fields
forAll(specieNames,i)
{
    //- Update mol concentration by mol fraction
    specieTable[specieNames[i]]->updateMolC();

    //- Initiate molPhi
    molPhi[i] = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "molPhi" + specieNames[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "molPhi",
                dimensionSet(0,0,-1,0,1,0,0),
                0.0
            )
        )
    );

    //- Initiate sourceTerm
    sourceTerm[i] = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sourceTerm" + specieNames[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "sourceTerm",
                dimensionSet(0,-3,-1,0,1,0,0),
                0.0
            )
        )
    );

    //- Initiate deltaC
    deltaC[i] = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "deltaC" + specieNames[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "deltaC",
                dimensionSet(0,-3,0,0,1,0,0),
                0.0
            )
        )
    );
}

//- Initiate phiConvg
phiConvg = - Mgf & (fvc::interpolate(fvc::grad(p) - rhog * g)) & mesh.Sf();


//- Relax factors
//- RelaxFactor for each gas species
List<scalar> relaxFactor(nSpecie);
scalar SwRelaxFactor = 1.0;
scalar TRelaxFactor = 1.0;

forAll(specieNames,i)
{
    relaxFactor[i] = 1.0;
}

//- Initialte thermalCapacity
thermalCapacitys =  (rhos + substrate) * thermalCaps;
thermalCapacityw = eps * Sw * rhow * thermalCapw;
thermalCapacityg = eps * (1-Sw) * rhog * thermalCapg;