label tmplduMatrix=lduMatrix::debug;
label tmpsolverPerf=solverPerformance::debug;
label tmpsolverPerfv=SolverPerformance<vector>::debug;
lduMatrix::debug = 0; 
solverPerformance::debug = 0;
SolverPerformance<vector>::debug = 0;
p = tp;
U = tU;
e = te;

solve
(
    fvm::ddt(p) == pFlux
);
p.correctBoundaryConditions();
dp = p - Pinf;

solve
(
    fvm::ddt(U) == UFlux
);
U.correctBoundaryConditions();

solve
(
    fvm::ddt(e) == thermo.Cv()*TFlux
);

e.correctBoundaryConditions();

lduMatrix::debug=tmplduMatrix;
solverPerformance::debug=tmpsolverPerf;
SolverPerformance<vector>::debug=tmpsolverPerfv;

thermo.correct();

rho.ref() = psi()*p.ref();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

rhoU.ref() = rho.ref()*U.ref();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

rhoE.ref() = p.ref()/(thermo.gamma()->ref()-1.0)
                             + 0.5*rho.ref()*magSqr(U.ref());
rhoE.boundaryFieldRef() ==
     p.boundaryField()/(thermo.gamma()->boundaryField()-1.0)
     + rho.boundaryField()*
     (
         0.5*magSqr(U.boundaryField())
     );

dT = T - Tinf;
