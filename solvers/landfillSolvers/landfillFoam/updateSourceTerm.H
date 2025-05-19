//- Update the source terms for the transport equations
methaneOxidation.update();
substrateDegrad->update();

forAll(specieNames,i)
{

    if(specieNames[i] == "O2")
    {
        sourceTerm[i].ref() = substrateDegrad->sourceO2() + methaneOxidation.sourceO2();
    }
    else if(specieNames[i] == "CH4")
    {
        sourceTerm[i].ref() =  substrateDegrad->sourceCH4() + methaneOxidation.sourceCH4();
    }
    else if(specieNames[i] == "CO2")
    {
        sourceTerm[i].ref() = substrateDegrad->sourceCO2() + methaneOxidation.sourceCO2();
    }
    else
    {
        sourceTerm[i].ref() = dimensionedScalar(dimMoles/dimVolume/dimTime, 0.0);
    }
}

sourceTermw =  (substrateDegrad->sourcew() + methaneOxidation.sourcew())/ rhow ;

heatSource = substrateDegrad->sourceHeat() + methaneOxidation.sourceHeat();