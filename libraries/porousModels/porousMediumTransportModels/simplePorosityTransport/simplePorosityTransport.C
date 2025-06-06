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

#include "fvCFD.H"
#include "simplePorosityTransport.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{
defineTypeNameAndDebug(simplePorosityTransport, 0);

addToRunTimeSelectionTable
(
    porousMediumTransportModel,
    simplePorosityTransport,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumTransportModels::simplePorosityTransport::simplePorosityTransport
(
    const word& phaseName,
    const fvMesh& mesh,
    const IOdictionary& transportProperties
)
    :
    porousMediumTransportModel(phaseName, mesh, transportProperties)
{}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::porousMediumTransportModels::simplePorosityTransport::solveTransport
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& theta,
    const volScalarField& exchangeTerm
)
{
    composition_.correct(U, theta);

    dictionary solverDict = composition_.Y(0).mesh().solver("C");
    forAll(composition_.Y(), speciesi)
    {
        auto& C = composition_.Y(speciesi);
        const auto& R = composition_.R(speciesi);
        const auto& Deff = composition_.Deff(speciesi);
        const auto& lambda = composition_.lambda(speciesi);
        const auto& sourceTerm_tracer = composition_.sourceTerm(speciesi);

        fvScalarMatrix CEqn
            (
                R * fvm::ddt(theta,C)
                + fvm::div(phi, C, "div(phi,C)")
                - fvm::laplacian(theta*Deff, C, "laplacian(Deff,C)")
                ==
                - sourceTerm_tracer
                - R * theta * fvm::Sp(lambda,C)
            );

        CEqn.solve(solverDict);
        info(speciesi, theta, phi);
    }
}

void Foam::porousMediumTransportModels::simplePorosityTransport::solveTransport
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& eps,
    const volScalarField& hwater,
    const volScalarField& seepageTerm,
    const scalar& zScale
)
{
    composition_.correct(U, eps);

    dictionary solverDict = composition_.Y(0).mesh().solver("C");
    forAll(composition_.Y(), speciesi)
    {
        auto& C = composition_.Y(speciesi);
        const auto& R = composition_.R(speciesi);
        const auto& Deff = composition_.Deff(speciesi);
        const auto& lambda = composition_.lambda(speciesi);
        const auto& sourceTerm = composition_.sourceTerm(speciesi);

        fvScalarMatrix CEqn
                (
                        eps* R * fvm::ddt(hwater,C)
                        + fvm::div(phi, C, "div(phi,C)")
                        - fvm::laplacian(eps * hwater * Deff, C, "laplacian(Deff,C)")
                        ==
                        - sourceTerm * dimensionedScalar(dimLength, zScale)
                        - eps * R * hwater * fvm::Sp(lambda,C)
                        - fvm::Sp(seepageTerm,C)
                );

        CEqn.solve(solverDict);
        volScalarField waterContent(eps*hwater);
        waterContent /= zScale;
        info(speciesi, waterContent, phi);
    }
}


void Foam::porousMediumTransportModels::simplePorosityTransport::info
(
    const label& speciei,
    const volScalarField& theta,
    const surfaceScalarField& phi
)
{
    const auto& C = composition_.Y(speciei);
    const auto& R = composition_.R(speciei);
    const auto& mesh = C.mesh();
    const auto& sourceTerm = composition_.sourceTerm(speciei);

    //- terminal display
    Info << "Concentration: Min(" << speciesNames_[speciei] << ") = " << gMin(C.internalField())
             << " Max(" << speciesNames_[speciei] << ") = " << gMax(C.internalField())
             << " mass(" << speciesNames_[speciei] << ") = " << fvc::domainIntegrate(R * C * theta).value()
             << endl;
    Info << speciesNames_[speciei] << " mass balance (kg/s): sourceTerm = " << fvc::domainIntegrate(sourceTerm).value() << " ; ";

    forAll(phi.boundaryField(),patchi)
    {
        if (mesh.boundaryMesh()[patchi].type() == "patch")
        {
            Info << phi.boundaryField()[patchi].patch().name() << " = " <<  gSum(phi.boundaryField()[patchi]*C.boundaryField()[patchi]) << " ; ";
        }
    }
    Info << endl;

}

// ************************************************************************* //
