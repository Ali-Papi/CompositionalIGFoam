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

#include "timestepManagerTruncation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //
Foam::timestepManagerTruncation::timestepManagerTruncation(
    const Time& runTime,
    const volScalarField& vf,
    const labelList* dryCells
)
    :
    runTime_(runTime),
    vf_(vf),
    dryCells_(dryCells),
    timeScheme_(vf.mesh().ddtScheme("ddt("+vf.name()+")")),
    d3dt3Operator_(vf.mesh(),runTime.deltaTValue()),
    d2dt2Operator_(vf.mesh()),
    dVmax_(0),
    dV2max_(0),
    dV3max_(0),
    Vmax_(0),
    V2max_(0),
    V3max_(0)
{
    //- Read truncation error value in controlDict
    truncationError_ = runTime.controlDict().getOrDefault<scalar>("truncationError_"+vf_.name(),
                                                                 runTime.controlDict().getOrDefault<scalar>("truncationError", GREAT));

    Info << nl << "Timestepping for field " << vf_.name() << + " is based on time-scheme truncation error with :"
    << nl << "{"
    << nl << "    truncationError = " << truncationError_
    << nl << "}"
    << endl;

    //- derivative initialization to keep 1st user-defined time step
    if (timeScheme_ == "backward")
    {
        V3max_ = max(SMALL, gMax(vf_.internalField()));
        dV3max_ = 3*truncationError_*(V3max_+VSMALL)/Foam::pow(runTime_.deltaTValue(),3);
    }
    else if  (timeScheme_ == "CrankNicolson")
    {
        V3max_ = max(SMALL, gMax(vf_.internalField()));
        dV3max_ = 12*truncationError_*(V3max_+VSMALL)/Foam::pow(runTime_.deltaTValue(),3);
    }
    else if (timeScheme_ == "Euler")
    {
        V2max_ = max(SMALL, gMax(vf_.internalField()));
        dV2max_ = 2*truncationError_*(V2max_+VSMALL)/Foam::pow(runTime_.deltaTValue(),2);
    }
    else if (timeScheme_ != "steadyState")
    {
        FatalErrorIn("timestepManagerTruncation.C") << "time scheme " << timeScheme_
            << " does not work with timestepManagerTruncation class" << exit(FatalError);
    }
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timestepManagerTruncation::~timestepManagerTruncation()
{}

// * * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * //

void timestepManagerTruncation::update1stOrder(bool dynamicMesh)
{
    volScalarField dVdT(vf_-vf_.oldTime());
    if (dryCells_) dryCellsDerivativeUpdate(dVdT);
    dVmax_ = 0;
    Vmax_ = 0;
    forAll(dVdT, celli)
    {
        if(mag(dVdT[celli]) > dVmax_)
        {
            Vmax_ = mag(vf_[celli]);
            dVmax_ = mag(dVdT[celli]);
        }
    }
    if (dynamicMesh) Vmax_ = gMax(vf_);
    else reduce(Vmax_, maxOp<scalar>());
    reduce(dVmax_, maxOp<scalar>());
    if (Vmax_ <= 0) Vmax_ = SMALL;
}

void timestepManagerTruncation::update2ndOrder(bool dynamicMesh)
{
    volScalarField dV2dT2(d2dt2Operator_.fvcD2dt2(vf_));
    if (dryCells_) dryCellsDerivativeUpdate(dV2dT2);
    dV2max_ = 0;
    V2max_ = 0;
    forAll(dV2dT2, celli)
    {
        if(mag(dV2dT2[celli]) > dV2max_)
        {
            V2max_ = mag(vf_[celli]);
            dV2max_ = mag(dV2dT2[celli]);
        }
    }
    if (dynamicMesh) V2max_ = gMax(vf_);
    else reduce(V2max_, maxOp<scalar>());
    reduce(dV2max_, maxOp<scalar>());
    if (V2max_ <= 0) V2max_ = SMALL;
}

void timestepManagerTruncation::update3rdOrder(bool dynamicMesh)
{
    volScalarField dV3dT3(d3dt3Operator_.fvcD3dt3(vf_));
    if (dryCells_) dryCellsDerivativeUpdate(dV3dT3);
    dV3max_ = 0;
    V3max_ = 0;
    forAll(dV3dT3, celli)
    {
        if(mag(dV3dT3[celli]) > dV3max_)
        {
            V3max_ = mag(vf_[celli]);
            dV3max_ = mag(dV3dT3[celli]);
        }
    }
    if (dynamicMesh) V3max_ = gMax(vf_);
    else reduce(V3max_, maxOp<scalar>());
    reduce(dV3max_, maxOp<scalar>());
    if (V3max_ <= 0) V3max_ = SMALL;
}

void timestepManagerTruncation::dryCellsDerivativeUpdate(volScalarField& field)
{
    const labelList& dryCells = *dryCells_;
    forAll(dryCells, pointi) field[dryCells[pointi]] = 0;
}

// * * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * //

void timestepManagerTruncation::updateDerivatives(bool dynamicMesh)
{
    update1stOrder(dynamicMesh);
    if (timeScheme_ == "Euler") update2ndOrder(dynamicMesh);
    else update3rdOrder(dynamicMesh);
}
  

scalar timestepManagerTruncation::computeTimestep
(
    const scalar& truncationError
)
{
    scalar dt = GREAT;
    if (timeScheme_ == "Euler")
    {
        dt =Foam::pow(2*truncationError*(V2max_+VSMALL)/(dV2max_+VSMALL),1./2.);
    }
    else if (timeScheme_ == "backward")
    {
        d3dt3Operator_.storeDeltaT00(runTime_.deltaT0Value());
        dt = Foam::pow(3*truncationError*(V3max_+VSMALL)/(dV3max_+VSMALL),1./3.);
    }
    else if (timeScheme_ == "CrankNicolson")
    {
        d3dt3Operator_.storeDeltaT00(runTime_.deltaT0Value());
        dt = Foam::pow(12*truncationError*(V3max_+VSMALL)/(dV3max_+VSMALL),1./3.);
    }
    reduce(dt, minOp<scalar>());
    return dt;
}

scalar timestepManagerTruncation::computeTimestep()
{
    scalar dt = computeTimestep(truncationError_);
    return dt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
