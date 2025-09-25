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


\*---------------------------------------------------------------------------*/
#include "OmegaTable.H"

const std::map<double,double> Foam::OmegaTable::OmegaTable={
                {	0.3	, 	2.649	}	,
                {	0.35	, 	2.468	}	,
                {	0.4	, 	2.314	}	,
                {	0.45	, 	2.182	}	,
                {	0.5	, 	2.066	}	,
                {	0.55	, 	1.965	}	,
                {	0.6	, 	1.877	}	,
                {	0.65	, 	1.799	}	,
                {	0.7	, 	1.729	}	,
                {	0.75	, 	1.667	}	,
                {	0.8	, 	1.612	}	,
                {	0.85	, 	1.562	}	,
                {	0.9	, 	1.517	}	,
                {	0.95	, 	1.477	}	,
                {	1	, 	1.44	}	,
                {	1.05	, 	1.406	}	,
                {	1.1	, 	1.375	}	,
                {	1.15	, 	1.347	}	,
                {	1.2	, 	1.32	}	,
                {	1.25	, 	1.296	}	,
                {	1.3	, 	1.274	}	,
                {	1.35	, 	1.253	}	,
                {	1.4	, 	1.234	}	,
                {	1.45	, 	1.216	}	,
                {	1.5	, 	1.199	}	,
                {	1.55	, 	1.183	}	,
                {	1.6	, 	1.168	}	,
                {	1.65	, 	1.154	}	,
                {	1.7	, 	1.141	}	,
                {	1.75	, 	1.128	}	,
                {	1.8	, 	1.117	}	,
                {	1.85	, 	1.105	}	,
                {	1.9	, 	1.095	}	,
                {	1.95	, 	1.085	}	,
                {	2	, 	1.075	}	,
                {	2.1	, 	1.058	}	,
                {	2.2	, 	1.042	}	,
                {	2.3	, 	1.027	}	,
                {	2.4	, 	1.013	}	,
                {	2.5	, 	1.0006	}	,
                {	2.6	, 	0.989	}	,
                {	2.7	, 	0.9782	}	,
                {	2.8	, 	0.9682	}	,
                {	2.9	, 	0.9588	}	,
                {	3	, 	0.95	}	,
                {	3.1	, 	0.9418	}	,
                {	3.2	, 	0.934	}	,
                {	3.3	, 	0.9267	}	,
                {	3.4	, 	0.9197	}	,
                {	3.5	, 	0.9131	}	,
                {	3.6	, 	0.9068	}	,
                {	3.7	, 	0.9008	}	,
                {	3.8	, 	0.8952	}	,
                {	3.9	, 	0.8897	}	,
                {	4	, 	0.8845	}	,
                {	4.1	, 	0.8796	}	,
                {	4.2	, 	0.8748	}	,
                {	4.3	, 	0.8703	}	,
                {	4.4	, 	0.8659	}	,
                {	4.5	, 	0.8617	}	,
                {	4.6	, 	0.8576	}	,
                {	4.7	, 	0.8537	}	,
                {	4.8	, 	0.8499	}	,
                {	4.9	, 	0.8463	}	,
                {	5	, 	0.8428	}	,
                {	6	, 	0.8129	}	,
                {	7	, 	0.7898	}	,
                {	8	, 	0.7711	}	,
                {	9	, 	0.7555	}	,
                {	10	, 	0.7422	}	,
                {	12	, 	0.7202	}	,
                {	14	, 	0.7025	}	,
                {	16	, 	0.6878	}	,
                {	18	, 	0.6751	}	,
                {	20	, 	0.664	}	,
                {	25	, 	0.6414	}	,
                {	30	, 	0.6235	}	,
                {	35	, 	0.6088	}	,
                {	40	, 	0.5964	}	,
                {	50	, 	0.5763	}	,
                {	75	, 	0.5415	}	,
                {	100	, 	0.518	}    };
