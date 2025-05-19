{
    Info << "Solving T equation" << endl;
    T.storePrevIter();
    deltaT.storePrevIter();


    fvScalarMatrix TEqn
    (
        //thermalCapacity * fvm::ddt(T)
        //+ T * fvc::ddt(thermalCapacity)
        fvm::ddt(thermalCapacitys, T)
        + fvm::ddt(thermalCapacityw, T)
        + fvm::ddt(thermalCapacityg, T)
        + rhow * thermalCapw * fvm::div(phiw, T)
        - fvm::laplacian(DT, T)
        ==
        heatSource
    );

    forAll(specieNames,i)
    {
        surfaceScalarField& molPhii = molPhi[i].ref();

        dimensionedScalar molWeighti = specieTable[specieNames[i]]->molWeight();
        dimensionedScalar cpmoli = specieTable[specieNames[i]]->cp() * molWeighti /1000;

        surfaceScalarField phiHeati_partialT = molPhii * cpmoli;
        fvScalarMatrix divHeatGas(fvm::div(phiHeati_partialT , T));
        TEqn += divHeatGas;
    } 



    TEqn.solve();

    deltaT = T - T.prevIter();


    if(iter == 1)
    {
        TRelaxFactor = 1.0;
    }
    else
    {
        TRelaxFactor *= autoExplicitRelaxFactor(T,deltaT);
    }

    Info << "RelaxFactor[T]:  " << TRelaxFactor << endl;

    T.relax(TRelaxFactor);
    T.correctBoundaryConditions();

    deltaT = T - T.prevIter();

    Info << "Maximum temperature = " << gMax(T.internalField()) << endl;
    Info << "Minimum temperature = " << gMin(T.internalField()) << endl;
}
