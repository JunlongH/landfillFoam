{
    Sw.storePrevIter();
    deltaSw.storePrevIter();

    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]) || isA<zeroGradientFvPatchField<vector> >(Uw.boundaryField()[patchi]))
        {
            phiwopc.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        
    }

    fvScalarMatrix SwEqn
    (
        fvm::ddt(eps, Sw)
        + fvc::div(phiwopc)
        + fvm::laplacian(dpcdSf * Mwf, Sw)
        == sourceTermw
    );

    SwEqn.solve();


    deltaSw = Sw - Sw.prevIter();

    {
        if(iter == 1)
        {
            SwRelaxFactor = 1.0;
        }
        else
        {
            SwRelaxFactor *= autoExplicitRelaxFactor(Sw,deltaSw);
        }

        Sw.relax(SwRelaxFactor);

        Info << "RelaxFactor[Sw]: " << SwRelaxFactor << endl;
    }

    Sw.correctBoundaryConditions();
    deltaSw = Sw - Sw.prevIter();

    //- Update water velocity
    phiw = SwEqn.flux() + phiwopc;

    Uw = fvc::reconstruct(phiw);

    Info << "Water Saturation: " << " Min(Sw) = " << gMin(Sw.internalField()) 
        << " Max(Sw) = " << gMax(Sw.internalField())  << endl;

}