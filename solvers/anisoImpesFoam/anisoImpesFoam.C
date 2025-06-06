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

Application
    anisoImpesFoam

Description
    Transient solver for incompressible two-phase flow (Darcy's law) in porous media
    using the IMPES method (IMplicit Pressure Explicit Saturation).
    Permeability is anisotropic (K == volTensorField)

\*---------------------------------------------------------------------------*/

#include "PMFversion.H"
#include "fvCFD.H"
#include "incompressiblePhase.H"
#include "capillarityModel.H"
#include "relativePermeabilityModel.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "patchEventFile.H"
#include "eventInfiltration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {  Foam::FatalError.exit(); }
    PMFversion solverV;

    Info << "Create time\n" << Foam::endl;
    Time runTime(Time::controlDictName, args);

    #include "createMesh.H"

    #include "createTimeControls.H"

    Info<< "\nReading g" << endl;
    const meshObjects::gravity& g = meshObjects::gravity::New(runTime);

    #include "createFields.H"
    #include "createSbFields.H"
    #include "readTimeControls.H"

    autoPtr<sourceEventFile> sourceEvent = sourceEventFile::New("sourceEventFileWater", transportProperties);
    sourceEvent->init(runTime, Sb.name(), mesh, sourceTerm.dimensions());
    autoPtr<outputEventFile> outputEvent = outputEventFile::New(runTime, mesh);
    outputEvent->addField(Sb, phi, eps, "waterMassBalance.csv", true);
    outputEvent->init();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (sourceEvent->isPresent()) sourceEvent->updateIndex(runTime.timeOutputValue());
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.timeOutputValue());
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);
        if (sourceEvent->isPresent())
        {
            sourceEvent->updateValue(runTime);
            sourceTerm = sourceEvent->dtValuesAsField();
        }

        //- Solve saturation equation (explicit)
        #include "SEqn.H"
        #include "updateSbProperties.H"

        //- Solve pressure equation (implicit)
        #include "pEqn.H"

        outputEvent->write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
