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

#include "porousMediumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::porousMediumTransportModel> Foam::porousMediumTransportModel::New
(
    const word& phaseName,
    const fvMesh& mesh,
    const IOdictionary& transportProperties
)
{
    word modelType(transportProperties.lookupOrDefault<word>("porousMediumModel", "simplePorosity"));
    modelType += "Transport";
    Info << nl << "Selecting porousTransportMedium model => " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
                "porousMediumTransportModel::New(...)"
            )   << "Unknown porousMediumTransportModel type "
                << modelType << nl << nl
                << "Valid porousMediumTransportModels are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }

    return autoPtr<porousMediumTransportModel>
        (cstrIter()(phaseName, mesh, transportProperties));
}


// ************************************************************************* //
