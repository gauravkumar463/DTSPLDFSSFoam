if(unsteady == true)
{
	dimensionedScalar tolp("rhoTolerance", dimless, 1e-3);
	residualConvergence.lookup(tolp.name()) >> tolp.value();

	dimensionedScalar tolU("rhoUTolerance", dimless, 1e-3);
	residualConvergence.lookup(tolU.name()) >> tolU.value();

	dimensionedScalar tole("rhoETolerance", dimless, 1e-3);
	residualConvergence.lookup(tole.name()) >> tole.value();

	bool checkAverageRes = false;
	residualConvergence.readIfPresent("checkAverageResiduals", checkAverageRes);

	bool noCheckU = false;
	bool noCheckV = false;
	bool noCheckW = false;
	residualConvergence.readIfPresent("DontCheckUx", noCheckU);
	residualConvergence.readIfPresent("DontCheckUy", noCheckV);
	residualConvergence.readIfPresent("DontCheckUz", noCheckW);

	if(checkAverageRes == false)
	{
		if(itr==0){
			if(tresp==0) tresp=SMALL;
			const Vector<label>& directions = mesh.geometricD();
			for (direction dir=0; dir<directions.nComponents; dir++)
			{
			    if(tresU.component(dir)==0) tresU.component(dir)= SMALL;
			}
			if(trese==0) trese=SMALL;
			iresp=tresp;
			iresU=tresU;
			irese=trese;
			Info << "Initial Residual rho = " << tresp << ", rhoU = " << tresU << ", rhoE = " << trese << endl;
			Info << "Initial Average residual rho = " << tresp1 << ", rhoU = " << tresU1 << ", rhoE = " << trese1 << endl;
		}

		if(
			  (
			   		((tresp/iresp <= tolp.value())) &&
					(!(tresU.component(0)/iresU.component(0) > tolU.value() && noCheckU==false)) && 
					(!(tresU.component(1)/iresU.component(1) > tolU.value() && noCheckV==false)) && 
					(!(tresU.component(2)/iresU.component(2) > tolU.value() && noCheckW==false)) && 
					((trese/irese <= tole.value()))
			  )
		  )
		{
			converged = true;
		}

		if(converged == true || itr == maxItr-1){
		    Info << "Final Residual rho = " << tresp/iresp << ", rhoU = " << cmptDivide(tresU,iresU) << ", rhoE = " << trese/irese << endl;
		    if(converged == true){
		    	Info << "Iteration converged in " << itr + 1 << " iterations" << endl;
		    	break;
			}
		}
		else
		{
			Info << "Residual rho = " << tresp/iresp << ", rhoU = " << cmptDivide(tresU,iresU) << ", rhoE = " << trese/irese << endl;
		} 
	}
	else
	{
		if(itr==0){
			if(tresp1==0) tresp1=SMALL;
			const Vector<label>& directions = mesh.geometricD();
			for (direction dir=0; dir<directions.nComponents; dir++)
			{
			    if(tresU1.component(dir)==0) tresU1.component(dir)= SMALL;
			}
			if(trese1==0) trese1=SMALL;
			iresp=tresp1;
			iresU=tresU1;
			irese=trese1;
			Info << "Initial Residual rho = " << tresp << ", rhoU = " << tresU << ", rhoE = " << trese << endl;
			Info << "Initial Average residual rho = " << tresp1 << ", rhoU = " << tresU1 << ", rhoE = " << trese1 << endl;
		}

		if(
			  (
			   		((tresp1/iresp <= tolp.value())) &&
					(!(tresU1.component(0)/iresU.component(0) > tolU.value() && noCheckU==false)) && 
					(!(tresU1.component(1)/iresU.component(1) > tolU.value() && noCheckV==false)) && 
					(!(tresU1.component(2)/iresU.component(2) > tolU.value() && noCheckW==false)) && 
					((trese1/irese <= tole.value()))
			  )
		  )
		{
			converged = true;
		}

		if(converged == true || itr == maxItr-1){
		    Info << "Final Residual rho = " << tresp1/iresp << ", rhoU = " << cmptDivide(tresU1,iresU) << ", rhoE = " << trese1/irese << endl;
		    if(converged == true){
		    	Info << "Iteration converged in " << itr + 1 << " iterations" << endl;
		    	break;
			}
			Info << endl;
		}
		else
		{
			Info << "Residual rho = " << tresp1/iresp << ", rhoU = " << cmptDivide(tresU1,iresU) << ", rhoE = " << trese1/irese << endl;
		} 
	}
}
else
{
	Info << "Residual rho = " << tresp << ", rhoU = " << tresU << ", rhoE = " << trese << endl;	
	Info << "Average residual rho = " << tresp1 << ", rhoU = " << tresU1 << ", rhoE = " << trese1 << endl;
}
