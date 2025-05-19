Info << nl << "Reading gasSpecieDict" << endl;

IOdictionary gasSpecieDict
(
    IOobject
    (
        "gasSpecieDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- Construction of gasCompositionSpace
gasCompositionSpace gasComp
(
    mesh,
    p,
    T,
    gasSpecieDict
);

//- Reference to gas composition space members
HashPtrTable<fluidGasSpecie,word>& specieTable = gasComp.gasSpecieTable();
const label nSpecie = gasComp.nGasSpecie();
const wordList& specieNames = gasComp.gasSpecieNames();

const volScalarField& rhog = gasComp.massConcentration();
const volScalarField& thermalCapg = gasComp.cp();

//const surfaceScalarField& phig = gasComp.phig();


//- List of deltaCi
List<tmp<volScalarField>> deltaC(nSpecie);

//- List of gas specie mole diffusion fluxes
List<tmp<surfaceScalarField>> molPhiDiff(nSpecie);

//- List of gas specie mole total fluxes
List<tmp<surfaceScalarField>> molPhi(nSpecie);

//- List of sourceTerm for gases and water
List<tmp<volScalarField>> sourceTerm(nSpecie);


Info << "Creating Dij" << endl;
//- Diffusion coefficients
Foam::porousBinaryDiffusionCoeffField Dij
(
    mesh,
    p,
    T,
    eps,
    Sg,
    gasComp
);

Info << "Creating DiK" << endl;
Foam::knudsenDiffusionCoeffField DiK
(
    mesh,
    T,
    kgmag,
    kmag,
    gasComp,
    "Heid"
);

