{
    volScalarField& rDeltaT = trSubDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    rDeltaT.ref().field() = sumAmaxSf/(nCoNum*mesh.V());
    if(preconditioning==false){
    	rDeltaT.ref() = 1.0/phyDeltaT;
    }
    else{
    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
