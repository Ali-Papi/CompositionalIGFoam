{
    
    forAll(Sb, cellI)
    {
	    if (Sb[cellI]>= Sb0[cellI]-1e-5 && Sb[cellI]<= Sb0[cellI]+1e-5) 
	    {
		    yCO2[cellI]= 0;
		    xCO2[cellI]= 0;
	    }
	    else
	    {
			    xCO2[cellI] = (p[cellI]-Psat.value())/(H.value()-Psat.value());
			    yCO2[cellI] = xCO2[cellI] / PH[cellI];
	    }
    }

    xH2O = 1-xCO2;
    yH2O = 1-yCO2;


//    Info << "xH2O: " << " Min(xH2O) = " << gMin(xH2O) << " Max(xH2O) = " << gMax(xH2O) << endl;

}
