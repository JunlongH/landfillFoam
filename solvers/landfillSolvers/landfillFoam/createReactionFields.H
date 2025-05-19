//////////////////////////////////////////////////////////////////
////////////////////// REACTION PROPERTIES ///////////////////////
//////////////////////////////////////////////////////////////////
Info << "Creating substrate concentration field " << endl;
volScalarField substrate
(
    IOobject
    (
        "substrate",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "substrate",
        dimDensity,
        50 //kg/m3
    )
);

Info << nl << "Reading reactionProperties" << endl;

IOdictionary reactionProperties
(
    IOobject
    (
        "reactionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<landfillBiodegradation> substrateDegrad = landfillBiodegradation::New
(
    "landfillBiodegradation",
    reactionProperties,
    substrate
);

methaneOxidation methaneOxidation
(
    reactionProperties,
    T,
    p,
    Sw
);


