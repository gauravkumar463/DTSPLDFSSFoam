bool unsteady = true;
scalar dtSwitch = 1.0;
mesh.schemesDict().readIfPresent("unsteady", unsteady);

bool inviscid = false;
if (mesh.schemesDict().readIfPresent("inviscid", inviscid))
{
    if (inviscid == true)
    {
        Info<< "Running solver in inviscid mode" << endl;
    }
}

bool preconditioning = true;
if (mesh.schemesDict().readIfPresent("preconditioning", preconditioning))
{
    if (preconditioning == true)
    {
	if(unsteady == true){
            Info<< "Running unsteady solver with preconditioning" << endl;
	}
	else{
            Info<< "Running steady solver with preconditioning" << endl;
	    dtSwitch = 0.0;
	}
    }
    else{
	dtSwitch = 0.0;
	if(unsteady == true){
            Info<< "Running unsteady solver without preconditioning" << endl;
	}
	else{
            Info<< "User input: steady solver. But, running unsteady solver without preconditioning" << endl;
	}
    }
}
