const dictionary& IterationDict = mesh.solution().subOrEmptyDict("Iteration");
const scalar CiResidualTol = IterationDict.lookupOrDefault<scalar>("tolerance", GREAT);
const scalar TResidualTol = IterationDict.lookupOrDefault<scalar>("Ttolerance", CiResidualTol);
const label maxIter = IterationDict.lookupOrDefault<label>("maxIter",10);

Info << nl << "Iteration loop control" << nl << "{"
<< nl << "    tolerance = " << CiResidualTol
<< "    maximum number of iteration = " << maxIter
<< nl << "}" << endl;


if (!IterationDict.found("tolerance") )
{
    FatalErrorIn("readIterationControls.H")
        << "Missing keyword tolerance in system/fvSolution.Iteration"
            << abort(FatalError);
}
