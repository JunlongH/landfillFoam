forAll(specieNames,i)
{
    volScalarField& Ci = specieTable[specieNames[i]]->molConcentration();

    if(gMin(Ci) < 0.0)
    {
        Ci = pos(Ci) * Ci;
    }
}

//- Restrain temperature to be above 298.15K
if(gMin(T.internalField()) < 298.15)
{
    /*T = pos(T-dimensionedScalar(dimTemperature,298.15)) * ( T - dimensionedScalar(dimTemperature,298.15) ) + dimensionedScalar(dimTemperature,298.15);
    volScalarField::Internal& Tref = T.ref();
    forAll(Tref, i)
    {
        if(Tref[i] < 298.15)
        {
            Tref[i] = 298.15;
        }
    }*/
   T = max(T,dimensionedScalar(dimTemperature,298.15));
    
}

gasComp.update();