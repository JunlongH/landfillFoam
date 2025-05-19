{
volScalarField& deltaCi =deltaC[i].ref();


deltaCi.storePrevIter();
deltaCi = Ci - Ci.prevIter();
if(iter == 1)
{
    relaxFactor[i] = 1.0;
}
else
{
    relaxFactor[i] *= autoExplicitRelaxFactor(Ci,deltaCi);
}


Ci.relax(relaxFactor[i]);
deltaCi = Ci - Ci.prevIter();

Info << "relaxFactor[" << specieNames[i] << "] = " << relaxFactor[i] << endl;

}