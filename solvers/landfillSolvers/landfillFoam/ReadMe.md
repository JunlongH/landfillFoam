# landfillAerationFoam_v1.2

Comparing with v1: 
1. Deleted water evaporation module.
2. Reordered field declaration and initialization
3. Modified convective heat transport using mole flux of each gas
```c++
    fvScalarMatrix TEqn
    (
        //thermalCapacity * fvm::ddt(T)
        //+ T * fvc::ddt(thermalCapacity)
        fvm::ddt(thermalCapacity, T)
        + rhow * thermalCapw * fvm::div(phiw, T)
        - fvm::laplacian(DT, T)
        ==
        heatSource
    );



    forAll(specieNames,i)
    {
        surfaceScalarField& molPhiDiffi = specieTable[specieNames[i]]->molPhiDiff();
        volScalarField& Ci = specieTable[specieNames[i]]->molConcentration();
        dimensionedScalar molWeighti = specieTable[specieNames[i]]->molWeight();
        dimensionedScalar cpmoli = specieTable[specieNames[i]]->cp() * molWeighti /1000;

        surfaceScalarField phiHeati_partialT = molPhiDiffi * cpmoli + phiConvg * fvc::interpolate(Ci) * cpmoli;
        fvScalarMatrix divHeatGas(fvm::div(phiHeati_partialT , T));
        TEqn += divHeatGas;
    } 
```

# landfillAerationFoam_v1.3

1. Addressed the problem that temperature be lower than 298.15
> This is because water flux and gas flux on faces be in consistent between mass transfer and heat transfer
> In mass transfer, water transport equation: 
```c++
    fvScalarMatrix SwEqn
    (
        fvm::ddt(eps, Sw)
        + fvc::div(phiwopc)
        + fvm::laplacian(dpcdSf * Mwf, Sw)
        == sourceTermw
    );
```
>Original water transport in heat transfer be:
```c++
phiGw = rhow * (Mwf & g) & mesh.Sf();
phiPw = (Mwf & gradpf) & mesh.Sf();
phiwopc = - phiPw + phiGw;
phipc = dpcdSf * (Mwf & fvc::interpolate(fvc::grad(Sw))) & mesh.Sf();
phiw = phiwopc + phipc;
```
> It's better to use: which is more consistent with SwEqn
```c++
phiw = SwEqn.flux() + phiwopc;
```
>  
> While for GAS transport
```c++
    fvScalarMatrix CiEqn
    (
        eps * fvm::ddt(Ci)
        - eps * fvm::ddt(Sw,Ci)
        + fvc::div(molPhiiDiff_)
        + fvm::div(phiConvg, Ci)
        == sourceTermi
    );
```
>Original:
```c++
    forAll(specieNames,i)
    {
        volScalarField& Ci = specieTable[specieNames[i]]->molConcentration();
        dimensionedScalar molWeighti = specieTable[specieNames[i]]->molWeight();
        dimensionedScalar cpmoli = specieTable[specieNames[i]]->cp() * molWeighti /1000;

        surfaceScalarField phiHeati_partialT = molPhiDiffi * cpmoli + phiConvg * fvc::interpolate(Ci) * cpmoli;
        fvScalarMatrix divHeatGas(fvm::div(phiHeati_partialT , T));
        TEqn += divHeatGas;
    } 
```
In this code, ```phiConvg * fvc::interpolate(Ci) * cpmoli``` might be in consistent with ```CiEqn```, depending on ```divSchemes```
>Should be:
```c++
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

    molPhi[i] = CiEqn.flux() + molPhiiDiff_;

    ......
    ......

    forAll(specieNames,i)
    {
        surfaceScalarField& molPhii = molPhi[i].ref();
        dimensionedScalar molWeighti = specieTable[specieNames[i]]->molWeight();
        dimensionedScalar cpmoli = specieTable[specieNames[i]]->cp() * molWeighti /1000;

        surfaceScalarField phiHeati_partialT = molPhii * cpmoli;
        fvScalarMatrix divHeatGas(fvm::div(phiHeati_partialT , T));
        TEqn += divHeatGas;
    } 

```

***Be aware that ```.flux()``` method can prevent surface flux be consistent whatever ```divSchemes``` be.***
