/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

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

Class
    Foam::sourceEventFile

Description
    Event file handler/reader for the porousMultiphaseFoam toolbox which
    contain 'm' dates with 'n' source points where flow rate is expressed 
    in m3/s. Source terms are applied in the closest cells to the given point
    The class reads the file :

    date time1
    x1 y1 z1 sourceTerm11
    x2 y2 z2 sourceTerm12
    ...
    xn yn zn sourceTerm1n
    date time2
    x1 y1 z1 sourceTerm21
    x2 y2 z2 sourceTerm22
    ...
    date timem
    ...
    xn yn zn sourceTermmn

SourceFiles
    sourceEventFile.C

\*---------------------------------------------------------------------------*/
#ifndef sourceEventFile_H
#define sourceEventFile_H

#include "fvMesh.H"
#include "eventFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluidsourceEventFile Declaration
\*---------------------------------------------------------------------------*/

class sourceEventFile
:
    public eventFile
{

private:

    label ncoordinates_;
    List<point> coordinates_;

    labelList idCoordinates_;
    dimensionSet fieldDims_ = dimless;

protected:
    void onMeshChanged() override;

public:

    // Constructors    
    //- Copy constructor
    sourceEventFile
    (
        const sourceEventFile& fileToCopy
    );

    //- Construct from Istream
    sourceEventFile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    ~sourceEventFile();

    // Member Functions
    
    //- Access function
    const label& ncoordinates() const {return ncoordinates_;}
    const List<point>& coordinates() const {return coordinates_;}


    //- Access as a volScalarField
    void setFieldDimensions(const dimensionSet& dims);
    tmp<volScalarField> dtValuesAsField() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
