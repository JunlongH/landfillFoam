/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is  an extension of OpenFOAM.

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

Application
    porousReactionHeatFoam

Description
    Transient solver for multi-component gas flow / two-phase flow/
    first order reaction and heat transfer .

\*---------------------------------------------------------------------------*/
//- OpenFoam original libraries
#include "fvMesh.H"
#include "Time.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "GeometricField.H"
#include "constants.H"
#include "argList.H"
#include "uniformDimensionedFields.H"
#include "fixedValueFvPatchField.H" 
#include "zeroGradientFvPatchField.H"

//- User defined libraries
#include "gasCompositionSpace.H" //gas composition
#include "fluidGasSpecie.H" //gas specie
#include "porousBinaryDiffusionCoeffField.H" //binary diffusion coefficient in porous media
#include "knudsenDiffusionCoeffField.H" //knudsen diffusion coefficient
#include "dustyGas.H" //dusty gas model
#include "maxwellStefan.H" //maxwell stefan model
#include "incompressiblePhase.H" //incompressible phase model for water
#include "relativePermeabilityModel.H" //relative permeability model
#include "capillarityModel.H" //capillarity model
#include "waterEvaporationModel.H" //water evaporation model
#include "autoExplicitRelaxFactor.H" //function for calculating relaxation factor for each field
#include "checkResidual.H"  //function for checking residuals for each field
#include "landfillBiodegradation.H" //landfill biodegradation model
#include "methaneOxidation.H" //methane oxidation model


#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createSwFields.H"
    #include "createThermalFields.H"
    #include "createSpecies.H"
    #include "createReactionFields.H"

    #include "initiateFields.H"
    #include "readIterationControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    label iter=0;

    while (runTime.run())
    {
        #include "modifyTimeStep.H"
        Info<< "Time = " << runTime.timeName() << nl << endl;

        iter = 0;
        bool converge = false;

        #include "updateProperties.H"
        
        while ( !converge && iter < maxIter )
        {
            iter++;
            Info << "Time: "<< runTime.timeName() << " Iteration " << iter << endl;
            
            #include "updateSourceTerm.H"
            #include "cEqn.H"
            #include "TEqn.H"
            #include "SwEqn.H"
            #include "updateProperties.H"
            #include "computeResiduals.H"
        }
        if (!converge && iter >= maxIter)
        {
            Info << endl;
            FatalErrorIn("landfillFoam.C") << "Non-convergence with fixed timestep => Decrease the time step or increase tolerance" 
            << exit(FatalError);
        }

        Info << "Pressure p: " << " Min(p) = " << gMin(p.internalField()) << " Max(p) = " << gMax(p.internalField()) << endl;

        #include "checkZero.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" 
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
