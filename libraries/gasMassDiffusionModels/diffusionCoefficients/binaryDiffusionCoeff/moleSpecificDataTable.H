/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Data
    Foam::moleSpecificDataTable

Description
    A map data storing the molecular specific data for computing binary diffusion
    coefficients.

SourceFiles
    molSpecificDataTable.C

\*---------------------------------------------------------------------------*/

#ifndef moleSpecificDataTable_H
#define moleSpecificDataTable_H

#include <array>
#include "HashTable.H"
#include "scalar.H"
#include "word.H"

namespace Foam
{
    namespace moleSpecificDataTable
    {
        const HashTable<std::array<scalar,3>, word> moleSpecificDataTable=
        {
            {"H2",{2.016,2.915,38.}},
            {"He",{4.003,2.576,10.2}},
            {"Ne",{20.183,2.789,35.7}},
            {"Ar",{39.948,3.432,122.4}},
            {"Kr",{83.80,3.675,170.0}},
            {"Xe",{131.30,4.009,234.7}},
            {"Air",{28.97,3.617,97}},
            {"CO2",{44.01,3.996,190.}},
            {"CO",{28.01,3.59,110.}},
            {"N2",{28.013,3.667,99.8}},
            {"H2O",{18.016,2.641,809.1}},
            //{"CH4",{16.042,3.758,148.6}},
            {"O2",{32.,3.433,113.}},
            {"NO",{30.01,3.470,119.0}},
            {"N2O",{44.01,3.879,220.0}},
            {"SO2",{64.06,4.026,363.0}},
            {"F2",{38.00,3.653,112.0}},
            {"Cl2",{70.91,4.115,357.0}},
            {"Br2",{159.82,4.268,520.0}},
            {"I2",{253.81,4.982,550.0}},
            {"CH4",{16.04,3.780,154.0}},
            //{"acetylene",26.04,4.114,212.0},
            //{"Subscript[CH, 2]==Subscript[CH, 2]",28.05,4.228,216.0},
            //{"Subscript[C, 2]Subscript[H, 6]",30.07,4.388,232.0},
            //{"Subscript[CH, 3]CH\[Congruent]CH",40.06,4.742,261.0},
            //{"Subscript[CH, 3]CH==Subscript[CH, 2]",42.08,4.766,275.},
            //{"Subscript[C, 3]Subscript[H, 6]",44.10,4.934,273.0},
            //{"n-Subscript[C, 4]Subscript[H, 10]",58.12,5.604,304.0},
            //{"i-Subscript[C, 4]Subscript[H, 10]",58.12,5.393,295.0},
            //{"n-Subscript[C, 5]Subscript[H, 12]",72.15,5.850,326.0},
            //{"i-Subscript[C, 5]Subscript[H, 12]",72.15,5.812,327.0},
            //{"C(Subscript[CH, 3]Subscript[), 4]",72.15,5.759,312.0},
            //{"n-Subscript[C, 6]Subscript[H, 14]", 86.18,6.264,342.0},
            //{"n-Subscript[C, 7]Subscript[H, 16]", 100.20,6.663,352.0},
            //{"n-Subscript[C, 8]Subscript[H, 18]", 114.23,7.035,361.0},
            //{"n-Subscript[C, 9]Subscript[H, 20]", 128.26,7.463,351.0},
            //{"Cyclohexane", 84.16,6.143,313.0},
            {"Benzene", {78.11,5.443,387.0}}
        };
    };

} // namespace Foam


#endif