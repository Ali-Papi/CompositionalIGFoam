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
    Foam::dualStaticMesh

Description
    A class for handling fields transfer between dual meshes (the secondary
    mesh is static).

SourceFiles
    dualStaticMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dualStaticMesh_H
#define dualStaticMesh_H

#include "multiMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class porousMediumModel Declaration
\*---------------------------------------------------------------------------*/

class dualStaticMesh
:
    public multiMesh
{
private:

    // Protected data
    //- Secondary mesh
    autoPtr<dynamicFvMesh> fineMeshPtr_;

    //- Time reference required for initial refinement
    Time& runTime_;

    //- Pointer to the cell mapping between two meshes
    labelList cellMapping_;
    labelListList boundaryMapping_;

    //- If mesh has already be refined
    bool refined_;

    //- List of tuple_field
    Tuple2<List<volScalarField*>, List<volScalarField*>> scalarFields_;
    Tuple2<List<volVectorField*>, List<volVectorField*>> vectorFields_;
    Tuple2<surfaceScalarField*, surfaceScalarField*> phiFields_;

    //- Initial refinement of the mesh and mapping computation
    void initialRefinement();

    //- Function to map fields between dual meshes
    template<class Type, template<class> class PatchField>
        void mapFieldCoarseToFine(
            GeometricField<Type, PatchField, volMesh>& field1,
            GeometricField<Type, PatchField, volMesh>& field2
            );

    //- Disallow copy construct
    dualStaticMesh(const dualStaticMesh&);

    //- Disallow default bitwise assignment
    void operator=(const dualStaticMesh&);

public:

    //- Runtime type information
    TypeName("dualStaticMesh");

    //- Construct from components
    dualStaticMesh(Time& runTime, dynamicFvMesh& mesh);

    //- Destructor
    ~dualStaticMesh() {}

    //- Return mesh
    dynamicFvMesh& fineMesh() {return fineMeshPtr_.ref();}

    //- Increment list of dual-mesh fields
    volScalarField& addField(volScalarField&);
    volVectorField& addField(volVectorField&);
    surfaceScalarField& addField(surfaceScalarField&);

    bool dynamic() { return false;}
    bool update();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
