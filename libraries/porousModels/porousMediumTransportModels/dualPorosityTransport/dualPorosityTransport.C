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
#include "dualPorosityTransport.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMediumTransportModels
{
defineTypeNameAndDebug(dualPorosityTransport, 0);

addToRunTimeSelectionTable
(
    porousMediumTransportModel,
    dualPorosityTransport,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediumTransportModels::dualPorosityTransport::dualPorosityTransport
(
    const word& phaseName,
    const fvMesh& mesh,
    const IOdictionary& transportProperties
)
    :
    porousMediumTransportModel(phaseName, mesh, transportProperties),
    dualPorosityTransportCoeffs_(transportProperties_.subDict("dualPorosityCoeffs")),
    epsMatrix_
    (
        IOobject
        (
            "epsMatrix",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        transportProperties_.getOrDefault<dimensionedScalar>("epsMatrix", dimensionedScalar(dimless, 2))
    ),
    matrixComposition_(
        transportProperties_,
        speciesNames("Matrix"),
        mesh,
        word::null,
        &sourceEventList_,
        "C",
        dimless,
        "Matrix"
    ),
    a_(dualPorosityTransportCoeffs_.get<dimensionedScalar>("a")),
    beta_(dualPorosityTransportCoeffs_.get<dimensionedScalar>("beta")),
    gammaW_(dualPorosityTransportCoeffs_.get<dimensionedScalar>("gammaW")),
    alphaS_(gammaW_*matrixComposition_.Dm(0)*beta_/(a_*a_)),
    exchangeTermFromFracture_
    (
        IOobject
        (
            "exchangeTerm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("",dimless/dimTime,0)
    ),
    exchangeTermFromMatrix_("exchangeTermFromMatrix", exchangeTermFromFracture_),
    UMatrix_(transportProperties_.db().lookupObject<volVectorField>("U"+phaseName_+"Matrix")),
    phiMatrix_(transportProperties_.db().lookupObject<surfaceScalarField>("phiMatrix")),
    thetaMatrix_(transportProperties_.db().lookupObject<volScalarField>(phaseName_+"Matrix"))
{
    matrixComposition_.initRetardCoef(epsMatrix_);
}

// * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * * //

void Foam::porousMediumTransportModels::dualPorosityTransport::solveTransport
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& theta,
    const volScalarField& exchangeTerm
)
{
    //- compute tracer exchange term from water transfer
    forAll(exchangeTerm, celli)
    {
        if (exchangeTerm[celli] > 0)
        {
            exchangeTermFromFracture_[celli] = exchangeTerm[celli] ;
            exchangeTermFromMatrix_[celli] = 0;
        }
        else
        {
            exchangeTermFromFracture_[celli] = 0;
            exchangeTermFromMatrix_[celli] = exchangeTerm[celli];
        }
    }

    //- fracture part
    composition_.correct(U, theta);

    dictionary solverDict = composition_.Y(0).mesh().solver("C");
    forAll(composition_.Y(), speciesi)
    {
        auto& C = composition_.Y(speciesi);
        const auto& Cmatrix = matrixComposition_.Y(speciesi);
        const auto& R = composition_.R(speciesi);
        const auto& Deff = composition_.Deff(speciesi);
        const auto& lambda = composition_.lambda(speciesi);
        const auto& sourceTerm_tracer = composition_.sourceTerm(speciesi);

        fvScalarMatrix CEqn
            (
                R * fvm::ddt(theta, C)
                + fvm::div(phi, C, "div(phi,C)")
                - fvm::laplacian(theta*Deff, C, "laplacian(Deff,C)")
                ==
                - sourceTerm_tracer
                - R * theta * fvm::Sp(lambda,C)
                - fvm::Sp(exchangeTermFromFracture_, C)
                - exchangeTermFromMatrix_ * Cmatrix
                + alphaS_ * thetaMatrix_ * Cmatrix
                - alphaS_ * fvm::Sp(thetaMatrix_,C)
            );

        CEqn.solve(solverDict);

        Info << "Concentration: Min(" << speciesNames_[speciesi] << ") = " << gMin(C.internalField())
             << " Max(" << speciesNames_[speciesi] << ") = " << gMax(C.internalField())
             << " mass(" << speciesNames_[speciesi] << ") = " << fvc::domainIntegrate(C * theta).value()
             << endl;
    }

    //- matrix part
    matrixComposition_.correct(UMatrix_, thetaMatrix_);

    forAll(matrixComposition_.Y(), speciesi)
    {
        auto& C = matrixComposition_.Y(speciesi);
        const auto& Cfracture = composition_.Y(speciesi);
        const auto& R = matrixComposition_.R(speciesi);
        const auto& Deff = matrixComposition_.Deff(speciesi);
        const auto& lambda = matrixComposition_.lambda(speciesi);
        const auto& sourceTerm_tracer = matrixComposition_.sourceTerm(speciesi);

        fvScalarMatrix CMatrixEqn
            (
                R * fvm::ddt(thetaMatrix_, C)
                + fvm::div(phiMatrix_, C, "div(phi,C)")
                -  fvm::laplacian(thetaMatrix_* Deff, C, "laplacian(Deff,C)")
                ==
                - sourceTerm_tracer
                - thetaMatrix_ * R * fvm::Sp(lambda,C)
                + exchangeTermFromFracture_ * Cfracture
                + exchangeTermFromMatrix_ * C
                + alphaS_ * thetaMatrix_ * (Cfracture - C)
            );

        CMatrixEqn.solve(C.mesh().solver("C"));

        Info<< "Concentration: Min(" << speciesNames_[speciesi] << ") = " << gMin(C.internalField())
            << " Max(" << speciesNames_[speciesi] << ") = " << gMax(C.internalField())
            << " mass(" << speciesNames_[speciesi] << ") = " << fvc::domainIntegrate(R*C*thetaMatrix_).value()
            << endl;
    }
}


void Foam::porousMediumTransportModels::dualPorosityTransport::solveTransport
        (
                const volVectorField& U,
                const surfaceScalarField& phi,
                const volScalarField& eps,
                const volScalarField& hwater,
                const volScalarField& seepageTerm,
                const scalar& zScale
        )
{
    FatalErrorIn("dualPorosityTransport::solveTransport")
        << "Dual porosity model is not available for 2D solvers" << abort(FatalError);
}

void Foam::porousMediumTransportModels::dualPorosityTransport::info
(
    const label& speciei,
    const volScalarField& theta,
    const surfaceScalarField& phi
)
{}
// ************************************************************************* //
