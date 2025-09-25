/*---------------------------------------------------------------------------*\
License
    This file is part of landfill modelling project, an extension of OpenFOAM
    modeling landfill aeration.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gasSpecie.H"
#include "atomicWeights.H"
#include "atomicDiffusionVolumes.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// - Construct from dictionary
Foam::gasSpecie::gasSpecie
(
    const fvMesh& mesh,
    const dictionary& gasSpecieDict,
    const word& gasSpecieName
)
    :
    mesh_(mesh),
    dict_(),
    name_(gasSpecieName),
    molWeight_(dimensionSet(1,0,0,0,-1,0,0),0.0),
    cp_(dimSpecificHeatCapacity,0.0)
{
    dict_ = gasSpecieDict.subDict("gasSpecie." + gasSpecieName);
    //- Directly read from dictionary or calculate from elemental composition from dictionary
    const word& inputType = dict_.lookup<word>("inputType");

    //- Directly read from dictionary
    if (inputType == "direct")
    {
        molWeight_ = dict_.lookupOrDefault<dimensionedScalar>("molWeight",dimensionedScalar("molWeight",dimensionSet(1,0,0,0,-1,0,0), 28.0));
        cp_ = dict_.lookupOrDefault<dimensionedScalar>("cp",dimensionedScalar("cp",dimSpecificHeatCapacity, 1040.0)); // J/kg/K
    }
    //- Calculate from elemental composition from dictionary
    else if (inputType == "composition")
    {
        //- Read elemental composition from dictionary using specieElement class
        List<specieElement> specieComposition = dict_.lookup("specieComposition");
        molWeight_.name() = "molWeight" + gasSpecieName;
        molWeight_.dimensions() = dimensionSet(1,0,0,0,-1,0,0);
        scalar molWeightVal = 0;

        forAll(specieComposition, i)
        {
            label nAtoms = specieComposition[i].nAtoms();
            const word& elementName = specieComposition[i].name();
            // molecular weight
            if (atomicWeights.found(elementName))
            {
                molWeightVal += atomicWeights[elementName] * nAtoms;
            }
            else
            {
                FatalErrorInFunction
                << "Element " << elementName << " not found in atomicWeights table"
                << nl << exit(FatalError);
            }

        }
        molWeight_.value() = molWeightVal;
        cp_ = dict_.lookupOrDefault<dimensionedScalar>("cp",dimensionedScalar("cp",dimSpecificHeatCapacity, 1040.0)); // J/kg/K
    }
    else
    {
        FatalErrorInFunction
        << "Unknown inputType " << inputType << "  for gas specie " << gasSpecieName
        << " in gasSpecieDict"
        << nl << exit(FatalError);
    }
}

// - Construct from specieComposition
Foam::gasSpecie::gasSpecie
(
    const fvMesh& mesh,
    const List<specieElement>& specieComposition,
    const word& gasSpecieName
)
    :
    mesh_(mesh),
    dict_(),
    name_(gasSpecieName)
{
    molWeight_.name() = "molWeight" + gasSpecieName;
    molWeight_.dimensions() = dimensionSet(1,0,0,0,-1,0,0);
    scalar molWeightVal = 0;
    forAll(specieComposition, i)
    {
        label nAtoms = specieComposition[i].nAtoms();
        const word& elementName = specieComposition[i].name();
        // molecular weight
        if (atomicWeights.found(elementName))
        {
            molWeightVal += atomicWeights[elementName] * nAtoms;
        }
        else
        {
            FatalErrorInFunction
            << "Element " << elementName << " not found in atomicWeights table"
            << nl << exit(FatalError);
        }
    }
    molWeight_.value() = molWeightVal;
    cp_ = dimensionedScalar("cp",dimSpecificHeatCapacity, 1000.0); //By default, J/kg/K 
}

// - Gas specie autoPtr
Foam::autoPtr<Foam::gasSpecie> Foam::gasSpecie::New
(
    const fvMesh& mesh,
    const dictionary& gasSpecieDict,
    const word& gasSpecieName
)
{
    return autoPtr<gasSpecie>
        (
            new gasSpecie(mesh, gasSpecieDict, gasSpecieName)
        );
}

// - Gas specie autoPtr
Foam::autoPtr<Foam::gasSpecie> Foam::gasSpecie::New
(
    const fvMesh& mesh,
    const List<specieElement>& specieComposition,
    const word& gasSpecieName
)
{
    return autoPtr<gasSpecie>
        (
            new gasSpecie(mesh, specieComposition, gasSpecieName)
        );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gasSpecie::~gasSpecie()
{}
