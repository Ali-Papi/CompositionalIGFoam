forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateValue(runTime);

//- compute source term
if (sourceEventIsPresent)
{
    sourceEvent.updateValue(runTime);
    sourceTerm = sourceEvent.dtValuesAsField();
}
//Info << "evt source term is = "<< sourceTerm << endl;