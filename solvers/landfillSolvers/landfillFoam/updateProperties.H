//- Gas saturation
Sg = 1.0 - Sw;

//- Update diffusion coefficients
Dij.update();
DiK.update();
//- Update gas diffusion flux
molPhiDiff =  gasMassDiffusion::dustyGas(mesh,gasComp,DiK.Dik(),Dij.DAB());


p.storePrevIter();

//- Update gas composition and pressure field within gasCompositionSpace class
gasComp.update();
p.correctBoundaryConditions();
deltap = p - p.prevIter();


//- Update gradp
gradpf = fvc::interpolate(fvc::grad(p));

//- Update permeability
krModel->correct();
kgmag = tr(K)/3 * krg;

krgf = fvc::interpolate(krg,"krg");
krwf = fvc::interpolate(krw,"krw");
Mgf = Kf*krgf/mug;
Mwf = Kf*krwf/muw;
Lwf = rhow*Kf*krwf/muw;


//- Update capillarity
pcModel->correct();
dpcdSf = fvc::interpolate(dpcdS,"dpcdS");

//- Update fluxes depending on saturation
phiGw = rhow * (Mwf & g) & mesh.Sf();
phiPw = (Mwf & gradpf) & mesh.Sf();
phiwopc = phiGw - phiPw;
//phipc = dpcdSf * (Mwf & fvc::interpolate(fvc::grad(Sw))) & mesh.Sf();
//phiw = phiwopc + phipc;
phiConvg = - Mgf & (gradpf - fvc::interpolate(rhog) * g) & mesh.Sf();

UConvG = fvc::reconstruct(phiConvg);


//- Update thermal capacity and diffusivity
thermalCapacitys =  (rhos + substrate + substrateDegrad->biomassA() + substrateDegrad->biomassN() + substrateDegrad->inertMass()) * thermalCaps; 
thermalCapacityw = eps * Sw * rhow * thermalCapw;
thermalCapacityg = dimensionedScalar(dimEnergy/dimTemperature/dimVolume,0.0);
forAll(specieNames,i)
{
    volScalarField& massCi = specieTable[specieNames[i]]->massConcentration();
    volScalarField thermalCapi = eps * Sg * specieTable[specieNames[i]]->cp()*massCi;
    thermalCapacityg += thermalCapi;
}


// average by volume
DT = DTs * (1 - eps) + DTw * eps * Sw + DTg * eps * Sg;