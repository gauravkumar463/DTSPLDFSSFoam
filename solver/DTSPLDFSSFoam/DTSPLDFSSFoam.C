/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "LDFSS/LDFSSFlux.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    // Courant numbers used to adjust the time-step
    scalarField sumAmaxSf;
    //label smoothItr=100;
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar iresp=0.0,irese=0.0;
    scalar tresp=0.0,trese=0.0,tresp1=0.0,trese1=0.0;
    vector iresU=vector::zero;
    vector tresU=vector::zero,tresU1=vector::zero;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        if(unsteady == true)
        {
            #include "updateDTSVariable.H"
        }

        #include "readTimeControls.H"
        #include "setDeltaT.H"
        
	residualConvergence = mesh.solutionDict().subDict("Convergence");
        scalar nCoNum = residualConvergence.lookupOrDefault<scalar>("CoNum",1); 

        dimensionedScalar phyDeltaT = runTime.deltaT();
        dimensionedScalar phyDeltaT0 = runTime.deltaT0();
	
        scalar maxItr = residualConvergence.lookupOrDefault<scalar>("maxItr",10);  
        
        TimeState subCycleTimeState = runTime.subCycle(1);
        label itr=0;
        bool converged = false;

        while(itr<maxItr)
        {
            //if((unsteady == true) && (preconditioning == true)){
              //  runTime++;
            //}
            Info << "\nIteration " << itr+1 << endl;
            
	    tp = p;
            tU = U;
            te = e;
            
	    for(label stage=0;stage<beta.size();stage++)
            {
                #include "LDFSS.H"
                #include "transformFluxes.H"

                if(stage == 0)
                {
		    if(preconditioning==true){
                    	#include "calculateResidual.H"
		    }
                    #include "centralCourantNo.H"
		    if(!(unsteady == true && preconditioning ==true)) nCoNum = CoNum;
                    #include "setRDeltaT.H"
                    trSubDeltaT.ref() = trSubDeltaT.ref()/beta[stage];
                }
                else
                {
                    trSubDeltaT.ref() = trSubDeltaT.ref()*beta[stage-1]/beta[stage];
                }
		
            	if((unsteady == true) && (preconditioning == true)){
                runTime++;
            	}
                #include "solveFluid.H"
            }
	    if(preconditioning==true){	
            	#include "checkConvergence.H"
	    }
            itr++;
        }
        runTime.endSubCycle();

        if(!inviscid == true){
            label nSubcycle = 1;
            TimeState subCycleTimeState = runTime.subCycle(nSubcycle);
            for(label i=0;i<nSubcycle;i++){
                turbulence->correct();
            }
            runTime.endSubCycle();
        }
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
