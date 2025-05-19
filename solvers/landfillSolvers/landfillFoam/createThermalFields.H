//////////////////////////////////////////////////////////////////
////////////////////// TEMPERATURE FIELD  ////////////////////////
//////////////////////////////////////////////////////////////////
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "T",
        dimTemperature,
        298.15
    )
);

volScalarField deltaT
(
    IOobject
    (
        "deltaT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);

volScalarField heatSource
(
    IOobject
    (
        "heatSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "heatSource",
        dimEnergy/dimTime/dimVolume,
        0.0
    )
);


Info<< "Reading heatTransferProperties\n" << endl;

IOdictionary heatTransferProperties
(
    IOobject
    (
        "heatTransferProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhos
(
    "rhos",
    dimDensity,
    heatTransferProperties.lookup("rhos")
);

dimensionedScalar thermalCaps
(
    "thermalCaps",
    dimSpecificHeatCapacity,
    heatTransferProperties.lookup("thermalCaps")
);

dimensionedScalar thermalCapw
(
    "thermalCapw",
    dimSpecificHeatCapacity,
    heatTransferProperties.lookup("thermalCapw")
);

dimensionedScalar DTs
(
    heatTransferProperties.lookupOrDefault("DTs",dimensionedScalar("DTs",dimEnergy/dimTemperature/dimTime/dimLength,3.990))
);

dimensionedScalar DTw
(
    heatTransferProperties.lookupOrDefault("DTw",dimensionedScalar("DTw",dimEnergy/dimTemperature/dimTime/dimLength,0.600))
);

dimensionedScalar DTg
(
    heatTransferProperties.lookupOrDefault("DTw",dimensionedScalar("DTw",dimEnergy/dimTemperature/dimTime/dimLength,0.025))
);


// Heat diffusion coefficient
volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "DT",
        dimEnergy/dimTemperature/dimTime/dimLength,
        0.095
    )
);



//- Total thermal capacity of solid phase
volScalarField thermalCapacitys
(
    IOobject
    (
        "thermalCapacity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimEnergy/dimTemperature/dimVolume
);

//- Total thermal capacity of liquid phase
volScalarField thermalCapacityw
(
    IOobject
    (
        "thermalCapacityw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimEnergy/dimTemperature/dimVolume
);

//- Total thermal capacity of gas phase
volScalarField thermalCapacityg
(
    IOobject
    (
        "thermalCapacityg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimEnergy/dimTemperature/dimVolume
);


