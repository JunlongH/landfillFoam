converge = true;
forAll(specieNames,i)
{

    volScalarField& Ci = gasComp.gasSpecieTable()[specieNames[i]]->molConcentration();

    volScalarField& deltaCi =deltaC[i].ref();

    converge = converge && checkResidual(Ci,deltaCi,CiResidualTol);
    Info << "Species " << specieNames[i] << " convergence: " << checkResidual(Ci,deltaCi,CiResidualTol) << endl;
}

converge = converge && checkResidual(Sw,deltaSw,TResidualTol)
    && checkResidual(p,deltap,TResidualTol)
    && checkResidual(T,deltaT,TResidualTol);

Info << "Sw convergence: " << checkResidual(Sw,deltaSw,TResidualTol) << endl;
Info << "p convergence: " << checkResidual(p,deltap,TResidualTol) << endl;
Info << "T convergence: " << checkResidual(T,deltaT,TResidualTol) << endl;
