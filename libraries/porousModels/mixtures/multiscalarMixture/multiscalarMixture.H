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
    Foam::multiscalarMixture

Description
    Multi-component mixture for porous transport.

SourceFiles
    multiscalarMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiscalarMixture_H
#define multiscalarMixture_H

#include "basicMultiComponentMixture.H"

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "typeInfo.H"

#include "dispersionModel.H"
#include "sourceEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiscalarMixture Declaration
\*---------------------------------------------------------------------------*/

class multiscalarMixture
:
    public basicMultiComponentMixture
{

public:

    //- Run time type information
    TypeName("multiscalarMixture");

    //- The base class of the mixture
    typedef multiscalarMixture basicMixtureType;

    // Constructors

    //- Construct from dictionary, species names, mesh and phase name
    multiscalarMixture
    (
        const dictionary& dict,
        const wordList& speciesNames,
        const fvMesh& mesh,
        const word& phaseName,
        List<sourceEventFile*>* sourceEventFileRegistry = nullptr,
        const word& sourceEventDtFieldNameOverride = "",
        const dimensionSet& sourceTermDimFactor = dimless,
        const word& porousRegion = ""
    );

    // Member Functions

    //- Update effective diffusivities and source terms
    void correct
    (
        const volVectorField& U,
        const volScalarField& theta
    );

    void correct
    (
        const volVectorField& U,
        const volScalarField& saturation,
        const volScalarField& eps
    );

    //- Retard coefficient [-]
    inline const volScalarField& R(const label speciesi) const;

    //- Decay coefficient [1/s]
    inline const dimensionedScalar& lambda(const label speciesi) const;

    //-  Molar diffusivity [m^2/s]
    inline const dimensionedScalar& Dm(const label speciesi) const;

    //- Effective diffusivity [m^2/s]
    inline const volTensorField& Deff(const label speciesi) const;

    //- Source term [Y(speciesi).dimensions()*sourceTermDimFactor/s]
    inline const volScalarField& sourceTerm(const label speciesi) const;

    //- List to the source event
    const PtrList<sourceEventFile>& sourceEvents() const;

    //- Source term [Y(speciesi).dimensions()*sourceTermDimFactor/s]
    inline const volScalarField& normalizedGradY() const;

    //- Compute retard coefficient and check eps
    bool initRetardCoef(const volScalarField& eps);

    //- Compute normalized maximal gradient
    void updateNormalizedGradY();

private:
    PtrList<volScalarField> R_;
    List<dimensionedScalar> lambdas_;
    PtrList<dispersionModel> dispersionModels_;
    PtrList<sourceEventFile> sourceEvents_;
    PtrList<volScalarField> sourceTerms_;
    List<dimensionedScalar> Kd_;
    List<dimensionedScalar> rs_;
    List<dimensionedScalar> epsTotal_;
    autoPtr<volScalarField> normalizedGradY_;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiscalarMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
