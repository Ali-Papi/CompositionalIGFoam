//- relative permeability (kr)
krModel->correct(Sb, true);
const volScalarField& kra = krModel->kra();
const volScalarField& krb = krModel->krb();
const volScalarField& dkradS = krModel->dkradS();
const volScalarField& dkrbdS = krModel->dkrbdS();

surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));
surfaceScalarField dkrafdS ("dkrafdS",fvc::interpolate(dkradS,"kra"));
surfaceScalarField dkrbfdS ("dkrbfdS",fvc::interpolate(dkrbdS,"krb"));

//- mobility and fractional flow 
surfaceScalarField Maf ("Maf",Kf*kraf/mua);
surfaceScalarField Laf ("Laf",rhoa*Kf*kraf/mua);	
surfaceScalarField Mbf ("Mbf",Kf*krbf/mub);
surfaceScalarField Lbf ("Lbf",rhob*Kf*krbf/mub);
surfaceScalarField Mf ("Mf",Maf+Mbf);
surfaceScalarField Lf ("Lf",Laf+Lbf);
surfaceScalarField Fbf ("Fbf",Mbf/Mf);
volScalarField Fb ("Fb",(krb/mub) / ( (kra/mua) + (krb/mub) ));

//- fluxes depending on saturation
pcModel->correct(Sb);
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
surfaceScalarField phiPc("phiPc",Mbf * fvc::interpolate(pcModel->dpcdS(),"pc")* fvc::snGrad(Sb) * mesh.magSf());
