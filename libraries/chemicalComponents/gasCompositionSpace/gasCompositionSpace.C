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

#include "gasCompositionSpace.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::gasCompositionSpace::gasCompositionSpace
(
    const fvMesh& mesh,
    volScalarField& p,
    const volScalarField& T,
    const dictionary& gasSpecieDict
)
    :
    name_(""),
    mesh_(mesh),
    p_(p),
    T_(T),
    C_
    (
        IOobject
        (
            "C"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionSet(0,-3,0,0,1,0,0)
    ),
    massC_
    (
        IOobject
        (
            "massC"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(1,-3,0,0,0,0,0)
    ),
    molWeight_
    (
        IOobject
        (
            "molWeight"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionSet(1,0,0,0,-1,0,0)
    ),
    cp_
    (
        IOobject
        (
            "thermalCapg"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimSpecificHeatCapacity
    ),
    phig_
    (
        IOobject
        (
            "phig"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phig", dimensionSet(0,3,-1,0,0,0,0), 0.0)
    ),
    /*U_
    (
        IOobject
        (
            "Ug",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionSet(0,1,-1,0,0,0,0)
    ),*/
    gasSpecieDict_(gasSpecieDict),
    gasSpecieNames_(gasSpecieDict_.lookup("gasSpecieNames"))
{
    // read gas species from dictionary
    
    nGasSpecie_ = gasSpecieNames_.size();
    forAll(gasSpecieNames_,i)
    {
        gasSpecieTable_.set
        (
            gasSpecieNames_[i],
            new fluidGasSpecie
            (
                mesh,
                gasSpecieDict_,
                gasSpecieNames_[i]+name_,
                p_,
                T_
            )
        );
    }

    //- initialize total mole concentration from gasSpecieTable_
    C_ = dimensionedScalar("C"+name_, dimensionSet(0,-3,0,0,1,0,0), 0.0);
    massC_ = dimensionedScalar("massC"+name_, dimensionSet(1,-3,0,0,0,0,0), 0.0);
    cp_ = dimensionedScalar("cp"+name_, dimSpecificHeatCapacity, 0.0);
    forAll(gasSpecieTable_,i)
    {
        C_ += gasSpecieTable_[gasSpecieNames_[i]]->molConcentration();
        massC_ += gasSpecieTable_[gasSpecieNames_[i]]->massConcentration();
    }
    
}

Foam::gasCompositionSpace::gasCompositionSpace
(
    word name,
    const fvMesh& mesh,
    volScalarField& p,
    const volScalarField& T,
    const dictionary& gasSpecieDict
)
    :
    name_(name),
    mesh_(mesh),
    p_(p),
    T_(T),
    C_
    (
        IOobject
        (
            "C"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionSet(0,-3,0,0,1,0,0)
    ),
    massC_
    (
        IOobject
        (
            "massC"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(1,-3,0,0,0,0,0)
    ),
    molWeight_
    (
        IOobject
        (
            "molWeight"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionSet(1,0,0,0,-1,0,0)
    ),
    cp_
    (
        IOobject
        (
            "thermalCapg"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimSpecificHeatCapacity
    ),
    phig_
    (
        IOobject
        (
            "phig"+name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phig", dimensionSet(0,3,-1,0,0,0,0), 0.0)
    ),
    /*U_
    (
        IOobject
        (
            "Ug",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionSet(0,1,-1,0,0,0,0)
    ),*/
    gasSpecieDict_(gasSpecieDict),
    gasSpecieNames_(gasSpecieDict_.lookup("gasSpecieNames"))
{
    // read gas species from dictionary
    
    nGasSpecie_ = gasSpecieNames_.size();
    forAll(gasSpecieNames_,i)
    {
        gasSpecieTable_.set
        (
            gasSpecieNames_[i],
            new fluidGasSpecie
            (
                mesh,
                gasSpecieDict_,
                gasSpecieNames_[i]+name_,
                p_,
                T_
            )
        );
    }

    //- initialize total mole concentration from gasSpecieTable_
    C_ = dimensionedScalar("C"+name_, dimensionSet(0,-3,0,0,1,0,0), 0.0);
    massC_ = dimensionedScalar("massC"+name_, dimensionSet(1,-3,0,0,0,0,0), 0.0);
    cp_ = dimensionedScalar("cp"+name_, dimSpecificHeatCapacity, 0.0);
    forAll(gasSpecieTable_,i)
    {
        C_ += gasSpecieTable_[gasSpecieNames_[i]]->molConcentration();
        massC_ += gasSpecieTable_[gasSpecieNames_[i]]->massConcentration();
    }
    
}

//- Destructor
Foam::gasCompositionSpace::~gasCompositionSpace()
{
    forAll(gasSpecieTable_,i)
    {
        delete gasSpecieTable_[gasSpecieNames_[i]];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

//- member functions
void Foam::gasCompositionSpace::updateMolConcentration()
{
    //- Update total mole concentration
    C_ = dimensionedScalar("C", dimensionSet(0,-3,0,0,1,0,0), 0.0);
    forAll(gasSpecieTable_,i)
    {
        C_ += gasSpecieTable_[gasSpecieNames_[i]]->molConcentration();
    }
} 

void Foam::gasCompositionSpace::updateMassConcentration()
{
    //- Update total mass concentration
    massC_ = dimensionedScalar("massC", dimensionSet(1,-3,0,0,0,0,0), 0.0);
    forAll(gasSpecieTable_,i)
    {
        massC_ += gasSpecieTable_[gasSpecieNames_[i]]->massConcentration();
    }
}

void Foam::gasCompositionSpace::updatePressure()
{
    //- Update pressure field
    p_ = C_ * Foam::constant::physicoChemical::R * T_;
}

void Foam::gasCompositionSpace::updateMolWeight()
{
    //- Update average molecular weight
    molWeight_ = dimensionedScalar(dimensionSet(1,0,0,0,-1,0,0), 0.0);
    forAll(gasSpecieTable_,i)
    {
        molWeight_ += gasSpecieTable_[gasSpecieNames_[i]]->volFrac() * gasSpecieTable_[gasSpecieNames_[i]]->molWeight();
    }
}

void Foam::gasCompositionSpace::updateCp()
{
    //- Update average specific heat capacity
    cp_ = dimensionedScalar("cp", dimSpecificHeatCapacity, 0.0);
    forAll(gasSpecieTable_,i)
    {
        cp_ += gasSpecieTable_[gasSpecieNames_[i]]->volFrac() * gasSpecieTable_[gasSpecieNames_[i]]->cp()* gasSpecieTable_[gasSpecieNames_[i]]->molWeight() 
            / (molWeight_+dimensionedScalar(dimensionSet(1,0,0,0,-1,0,0),VSMALL));  
    }
}

void Foam::gasCompositionSpace::updateVelocity()
{
    //- Update velocity field
    //U_ = dimensionedVector(dimensionSet(0,1,-1,0,0,0,0), Zero);
    //- Initialize average molecular weight
    phig_ = p_.db().lookupObject<surfaceScalarField>("phiConvg");
    forAll(gasSpecieTable_,i)
    {
    //    U_ += gasSpecieTable_[gasSpecieNames_[i]]->volFrac() * gasSpecieTable_[gasSpecieNames_[i]]->uSup() * gasSpecieTable_[gasSpecieNames_[i]]->molWeight() 
    //        / (molWeight_+dimensionedScalar(dimensionSet(1,0,0,0,-1,0,0),VSMALL));
        phig_ += gasSpecieTable_[gasSpecieNames_[i]]->molPhiDiff()*gasSpecieTable_[gasSpecieNames_[i]]->molWeight()*1e-3/(fvc::interpolate(massC_)+dimensionedScalar(dimensionSet(1,-3,0,0,0,0,0),VSMALL)) ;
    }

}

void Foam::gasCompositionSpace::update()
{

    //- Update total mole concentration
    updateMolConcentration();

    //- Update pressure field
    updatePressure();
    
    //- Note: update pressure should be done before update volumetric fraction
    //- Update all other fields for each gas specie
    forAll(gasSpecieTable_,i)
    {
        Info<< "Updating gas specie: " << gasSpecieNames_[i] << endl;
        gasSpecieTable_[gasSpecieNames_[i]]->updatebyMolC();
    }

    updateMassConcentration();

    //- Update molWeight should be before velocity
    //- Update average molecular weight
    updateMolWeight();

    //- Update average specific heat capacity
    //updateCp();
    
    //- Update velocity field
    //updateVelocity();

}