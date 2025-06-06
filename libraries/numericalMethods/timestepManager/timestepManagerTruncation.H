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
    Foam::fv::timestepManagerTruncation

Description
    Timestep manager compute timestep based on temporal truncation error. 
    Current version works with volScalarField and 3 time schemes:
    Euler / backward / CrankNicolson

SourceFiles
    timestepManagerTruncation.C

\*---------------------------------------------------------------------------*/

#ifndef timestepManagerTruncation_H
#define timestepManagerTruncation_H

#include "fvMesh.H"
#include "EulerD2dt2Scheme.H"
#include "EulerD3dt3Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class timestepManagerTruncation Declaration
\*---------------------------------------------------------------------------*/

class timestepManagerTruncation
{
    
private:

    //- reference to time
    const Time& runTime_;

    //- reference to field
    const volScalarField& vf_;

    //- User defined truncation error
    scalar truncationError_;

    //- reference to dry cells list
    const labelList* dryCells_;

    //- time scheme for this field
    const word timeScheme_;

    //- derivative operator
    EulerD3dt3Scheme d3dt3Operator_;
    fv::EulerD2dt2Scheme<scalar> d2dt2Operator_;

    //- maximal derivative/value
    scalar dVmax_, dV2max_, dV3max_;
    scalar Vmax_, V2max_, V3max_;

    void update1stOrder(bool);
    void update2ndOrder(bool);
    void update3rdOrder(bool);
    void dryCellsDerivativeUpdate(volScalarField&);
 
public:

    // Constructors

    //- Construct from field
    timestepManagerTruncation(
        const Time&,
        const volScalarField&,
        const labelList* dryCells
    );

    virtual ~timestepManagerTruncation();

    //- Access functions
    scalar Vmax(){return Vmax_;}
    scalar dVmax(){return dVmax_;}
    scalar V2max(){return V2max_;}
    scalar dV2max(){return dV2max_;}
    scalar V3max(){return V3max_;}
    scalar dV3max(){return dV3max_;}

    //- Access functions
    void updateDerivatives(bool dynamicMesh = false);
    scalar computeTimestep(const scalar&);
    scalar computeTimestep();
    

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
