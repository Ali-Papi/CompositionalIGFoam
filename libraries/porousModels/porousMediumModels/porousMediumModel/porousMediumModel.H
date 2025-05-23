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

Namespace
    Foam::porousMediumModel

Description
    A namespace for various relative permeability implementations.

Class
    Foam::porousMediumModel

Description
    An base class for relative permeability models.

SourceFiles
    porousMediumModel.C
    porousMediumModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef porousMediumModel_H
#define porousMediumModel_H

#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class porousMediumModel Declaration
\*---------------------------------------------------------------------------*/

class porousMediumModel
{
protected:

    // Protected data
    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Reference to the main dictionary
    const IOdictionary& transportProperties_;

    //- porosity field
    volScalarField eps_;

    //- permeability field
    volScalarField K_;

    //- injection/extraction source term
    volScalarField sourceTerm_;

    //- injection/extraction source term
    volScalarField exchangeTerm_;

    //- Disallow copy construct
    porousMediumModel(const porousMediumModel&);

    //- Disallow default bitwise assignment
    void operator=(const porousMediumModel&);

public:

    //- Runtime type information
    TypeName("porousMediumModel");

    // Constructors

    //- Construct from components
    porousMediumModel
    (
        const fvMesh& mesh,
        const IOdictionary& transportProperties,
        const dimensionSet& sourceDims,
        const word porousRegion = ""
    );

    //- Destructor
    virtual ~porousMediumModel()
        {}

    //- Return mesh
    const fvMesh& mesh() const {return mesh_;}

    //- Return transport properties
    const IOdictionary& transportProperties() const {return transportProperties_;}

    //- Return porosity
    const volScalarField& eps() const {return eps_;}

    //- Return permeability
    const volScalarField& K() const {return K_;}

    //- Return source term
    volScalarField& sourceTerm() {return sourceTerm_;}

    //- Return exchange term
    const volScalarField& exchangeTerm() const {return exchangeTerm_;}

    //- Check functions
    void check_eps() const;
    void check_K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
