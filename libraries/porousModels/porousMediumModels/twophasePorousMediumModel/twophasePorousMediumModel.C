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

\*---------------------------------------------------------------------------*/

#include "twophasePorousMediumModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(twophasePorousMediumModel, 0);
defineRunTimeSelectionTable(twophasePorousMediumModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twophasePorousMediumModel::twophasePorousMediumModel
(
    const word Sname,
    const fvMesh& mesh,
    const IOdictionary& transportProperties,
    const autoPtr<incompressiblePhase>& phase,
    const word porousRegion
)
    :
    porousMediumModel(mesh, transportProperties, dimless/dimTime, porousRegion),
    Sname_(Sname)
{
    pcModel_ = capillarityModel::New(mesh, transportProperties, Sname, porousRegion);
    krModel_ = relativePermeabilityModel::New(mesh, transportProperties, Sname, porousRegion);
}

// ************************************************************************* //
