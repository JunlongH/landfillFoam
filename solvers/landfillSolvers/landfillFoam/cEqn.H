
forAll(specieNames,i)
{
    volScalarField& Ci = specieTable[specieNames[i]]->molConcentration();
    Ci.storePrevIter();


    surfaceScalarField& molPhiiDiff_ = specieTable[specieNames[i]]->molPhiDiff();
    molPhiiDiff_ = molPhiDiff[i].ref();

    volScalarField& sourceTermi = sourceTerm[i].ref();

    fvScalarMatrix CiEqn
    (
        eps * fvm::ddt(Ci)
        - eps * fvm::ddt(Sw,Ci)
        + fvc::div(molPhiiDiff_)
        + fvm::div(phiConvg, Ci)
        == sourceTermi
    );

    CiEqn.solve();

    #include "relaxFactor.H"
    Ci.correctBoundaryConditions();

    //molPhi[i] = CiEqn.flux() + molPhiiDiff_;

    Info << "Max of " << specieNames[i] << " = " << gMax(Ci.internalField()) << endl;
    Info << "Min of " << specieNames[i] << " = " << gMin(Ci.internalField()) << endl;

}
