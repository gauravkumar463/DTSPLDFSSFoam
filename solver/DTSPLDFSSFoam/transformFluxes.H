rhoFlux = - fvc::div(phi);

rhoUFlux = - fvc::div(phiUp);

rhoEFlux = - fvc::div(phiEp);

if(!inviscid == true){
      volScalarField muEff("muEff", turbulence->muEff());
      volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));
      rhoUFlux += fvc::laplacian(muEff, U) + fvc::div(tauMC);// + fvOptions(rho,U);
      rhoEFlux += fvc::laplacian(turbulence->kappaEff(), T)
                      + fvc::div((
                            fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
                            + fvc::dotInterpolate(mesh.Sf(), tauMC)
                        ) & linearInterpolate(U));
}

if((unsteady == true) && (preconditioning == true))
{
    scalar coefft   = 1.0 + (phyDeltaT/(phyDeltaT + phyDeltaT0)).value();
    scalar coefft00 = (phyDeltaT*phyDeltaT/(phyDeltaT0*(phyDeltaT + phyDeltaT0))).value();
    scalar coefft0 = coefft + coefft00;

    if (mesh.moving())
    {
      rhoFlux.ref() -= ((coefft*V.internalField()*rho.internalField())-(coefft0*V0.internalField()*rho1.internalField())+(coefft00*V00.internalField()*rho2.internalField()))/(phyDeltaT*V.internalField());

      ddtRhoU.ref() = ((coefft*V.internalField()*rhoU.internalField())-(coefft0*V0.internalField()*rhoU1.internalField())+(coefft00*V00.internalField()*rhoU2.internalField()))/(phyDeltaT*V.internalField());
      ddtRhoU.correctBoundaryConditions();
      rhoUFlux.ref() -= ddtRhoU.internalField();

      rhoEFlux.ref() -= ((coefft*V.internalField()*rhoE.internalField())-(coefft0*V0.internalField()*rhoE1.internalField())+(coefft00*V00.internalField()*rhoE2.internalField()))/(phyDeltaT*V.internalField());
    }
    else
    {
      rhoFlux -= ((coefft*rho)-(coefft0*rho1)+(coefft00*rho2))/phyDeltaT;
      rhoUFlux -= ((coefft*rhoU)-(coefft0*rhoU1)+(coefft00*rhoU2))/phyDeltaT;
      rhoEFlux -= ((coefft*rhoE)-(coefft0*rhoE1)+(coefft00*rhoE2))/phyDeltaT;
    }
}

volScalarField delta = 1.5*dtSwitch/(phyDeltaT*trSubDeltaT.ref());

volScalarField theta =scalar(1.0)/Ur2+(thermo.gamma()-1.0)/c2;

volScalarField a1 = ((c2/thermo.gamma())*(1.0+theta*((U2/2.0)-(c2/(thermo.gamma()-1.0))))) 
                        + delta*((U2/2.0)-(c2/(thermo.gamma()*(thermo.gamma()-1.0))));

volScalarField a2 = scalar(1.0)/((theta*(c2/(thermo.gamma()-1.0)) -1.0) 
                                    + delta/(thermo.gamma()-1.0));

volScalarField a3 = (theta*c2/thermo.gamma()) + delta;

volScalarField eta = (1.0+delta)*(rho*thermo.Cv())
                            *(((thermo.gamma()-1.0)-(theta*c2))-delta);

pFlux = a2*((U2*rhoFlux/2.0) - (U & rhoUFlux) + rhoEFlux);

UFlux = (-(rhoFlux*U) + rhoUFlux)/(rho*(1.0+delta));

TFlux = (-a1*rhoFlux + a3*(U & rhoUFlux) - a3*rhoEFlux)/eta;
