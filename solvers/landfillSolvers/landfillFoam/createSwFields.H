//- relative permeability (kr)
krModel->correct();
const volScalarField& krg = krModel->kra();
const volScalarField& krw = krModel->krb();

volScalarField kgmag("kgmag",tr(K)/3*krg);


surfaceScalarField krgf ("krgf",fvc::interpolate(krg,"krg"));
surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw")); 

//- mobility computation 
surfaceTensorField Mgf ("Mgf",Kf*krgf/mug);
surfaceTensorField Mwf ("Mbf",Kf*krwf/muw);
surfaceTensorField Lwf ("Lbf",rhow*Kf*krwf/muw);

//- capillarity computation (pc)
pcModel->correct();
const volScalarField& dpcdS = pcModel->dpcdS();
surfaceScalarField dpcdSf("dpcdSf",fvc::interpolate(dpcdS,"dpcdS"));

//- fluxes depending on saturation
surfaceScalarField phiGw("phiGw",rhow * (Mwf & g) & mesh.Sf());

//- Water flux caused by p
surfaceScalarField phiPw("phiPw",(Mwf & gradpf) & mesh.Sf());

//- Water flux despite pc
surfaceScalarField phiwopc("phiwopc",- phiPw + phiGw);
surfaceScalarField phipc("phipc",dpcdSf * (Mwf & fvc::interpolate(fvc::grad(Sw))) & mesh.Sf());

//- Total water flux
surfaceScalarField& phiw = phasew->phi();
phiw = phiwopc + phipc;




