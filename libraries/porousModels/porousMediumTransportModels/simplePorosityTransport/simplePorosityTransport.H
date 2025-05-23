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
    Foam::porousMediumTransportModels::simplePorosityTransport

Class
    Foam::simplePorosityTransport

Description
    Simple porosity model for solute transport

SourceFiles
    simplePorosityTransport.C

\*---------------------------------------------------------------------------*/

#ifndef simplePorosityTransport_H
#define simplePorosityTransport_H

#include "volFields.H"
#include "porousMediumTransportModel.H"
#include "relativePermeabilityModel.H"
#include "capillarityModel.H"
#include "incompressiblePhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{

/*---------------------------------------------------------------------------* \
             Class simplePorosityTransport Declaration
\*---------------------------------------------------------------------------*/

class simplePorosityTransport
:
    public porousMediumTransportModel
{
protected:

    //- Disallow copy constructn
    simplePorosityTransport(const simplePorosityTransport&);

    //- Disallow default bitwise assignment
    void operator=(const simplePorosityTransport&);

public:

    //- Runtime type information
    TypeName("simplePorosityTransport");

    // Constructors

    //- Construct from components
    simplePorosityTransport
    (
        const word& phaseName,
        const fvMesh& mesh,
        const IOdictionary& transportProperties
    );

    //- Destructor
    ~simplePorosityTransport()
        {}

    //- Solve tracer transport equation (Richards' flow modeling)
    void solveTransport
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& theta,
        const volScalarField& exchangeTerm
    );

    //- Solve tracer transport equation (Dupuit-Darcy flow modeling)
    void solveTransport
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& eps,
        const volScalarField& hwater,
        const volScalarField& seepageTerm,
        const scalar& zScale
    );

    //- Display information to the terminal
    void info
    (
        const label& speciei,
        const volScalarField& theta,
        const surfaceScalarField& phi
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMediumTransportModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
