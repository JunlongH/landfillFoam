//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceVectorField gradpf("gradp", fvc::interpolate(fvc::grad(p),"gradp"));

Info << "Creating deltap field (for Iteration loop)" << endl;
volScalarField deltap
(
    IOobject
    (
        "deltap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p
);

Info << nl << "Reading field Sw" << endl;
volScalarField Sw
(
    IOobject
    (
        "Sw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Sg("Sg",1.0-Sw);

volScalarField deltaSw
(
    IOobject
    (
        "deltaSw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Sw
);



//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Swmin(transportProperties.lookupOrDefault("Swmin",dimensionedScalar("Swmin",dimless,0)));
dimensionedScalar Swmax(transportProperties.lookupOrDefault("Swmax",dimensionedScalar("Swmax",dimless,1)));

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//- Create water phase
autoPtr<incompressiblePhase> phasew = incompressiblePhase::New(mesh,transportProperties,"w");
const dimensionedScalar& rhow = phasew->rho();
const dimensionedScalar& muw = phasew->mu();
volVectorField& Uw = phasew->U();

dimensionedScalar mug
(
    transportProperties.lookupOrDefault
    (
        "mug",
        dimensionedScalar
        (
            "mug",
            dimDynamicViscosity,
            1.81e-5
        )
    )
);

//- Convective velocity
volVectorField UConvG
(
    IOobject
    (
        "UConvG",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(0,1,-1,0,0,0,0),vector::zero)
);


// solid phase rhog created in readHeatTransferProperties.H

volScalarField sourceTermw
(
    IOobject
    (
        "sourceTermw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sourceTermw", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

// Gas phase created in createSpecies.H


// Relative permeability 
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,Sw);
//scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));

// Capillary pressure models
autoPtr<capillarityModel> pcModel = capillarityModel::New("pcModel",transportProperties,Sw);

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Porosity
Info << nl << "Reading porosity field eps (if present)" << endl;
dimensionedScalar epsScalar(transportProperties.lookupOrDefault("eps",dimensionedScalar("",dimless,0.6)));
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    epsScalar
);

// Intrinsic permeability       
Info << nl << "Reading permeability field K" << endl;
volTensorField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// interpolated permeability
surfaceTensorField Kf(fvc::interpolate(K,"K"));

// magnitude of Intrinsic permeability
volScalarField kmag("kmag",tr(K)/3);

//- Gas flux caused by convection
surfaceScalarField phiConvg
(
    IOobject
    (
        "phiConvg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,3,-1,0,0,0,0)
);

