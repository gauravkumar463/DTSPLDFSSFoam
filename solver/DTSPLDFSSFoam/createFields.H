bool LTS = true;
#include "createRSubDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& p = thermo.p(); 

volScalarField& e = thermo.he();

volScalarField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("V", dimVolume, 1.0)
);
V.ref() = mesh.V();
V.storeOldTime();
volScalarField V0 = V.oldTime();
volScalarField V00 = V0.oldTime();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p/(thermo.gamma()-1.0) + 0.5*rho*magSqr(U)
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);

surfaceVectorField phiUp
(
    IOobject
    (
        "phiUp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (phi*linearInterpolate(U))
);

surfaceScalarField phiEp
(
    IOobject
    (
        "phiEp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi*linearInterpolate(e+0.5*magSqr(U))
);

surfaceScalarField amaxSf
(
    IOobject
    (
        "amaxSf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

volScalarField Ur2
(
    IOobject
    (
        "Ur2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
mesh,
dimensionedScalar("Ur2", dimVelocity*dimVelocity, 1.0),
"zeroGradient"
);

rho.storeOldTime();
volScalarField rho1 = rho.oldTime();
volScalarField rho2 = rho1.oldTime();

rhoU.storeOldTime();
volVectorField rhoU1 = rhoU.oldTime();
volVectorField rhoU2 = rhoU1.oldTime();

rhoE.storeOldTime();
volScalarField rhoE1 = rhoE.oldTime();
volScalarField rhoE2 = rhoE1.oldTime();

volScalarField rhoFlux
(
    IOobject
    (
        "rhoFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoFlux", dimensionSet(1,-3,-1,0,0,0,0), 1.0)
);

volVectorField rhoUFlux
(
    IOobject
    (
        "rhoUFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("rhoUFlux", dimensionSet(1,-2,-2,0,0,0,0), vector::one)
);

volScalarField rhoEFlux
(
    IOobject
    (
        "rhoEFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoEFlux", dimensionSet(1,-1,-3,0,0,0,0), 1.0)
);

volScalarField pFlux
(
    IOobject
    (
        "pFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pFlux", dimensionSet(1,-1,-3,0,0,0,0), 1.0)
);

volVectorField UFlux
(
    IOobject
    (
        "UFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("UFlux", dimensionSet(0,1,-2,0,0,0,0), vector::one)
);

volScalarField TFlux
(
    IOobject
    (
        "TFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TFlux", dimensionSet(0,0,-1,1,0,0,0), 1.0)
);

dictionary residualConvergence = mesh.solutionDict().subDict("Convergence");

//scalarList beta = {0.11, 0.2766, 0.5, 1.0};

label numberRungeKuttaStages = 1;
scalarList beta(numberRungeKuttaStages, 1.0);

if(mesh.solutionDict().found("multiStage"))
{
    // Read multi-stage coeffs dict
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    // Read number of multi-stages
    if (multiStage.found("numberRungeKuttaStages"))
    {
         multiStage.lookup("numberRungeKuttaStages") >> numberRungeKuttaStages;
         beta.setSize(numberRungeKuttaStages, 1.0);
    }

    // Read multi-stage coeffs
    if(multiStage.found("RKCoeff") && numberRungeKuttaStages > 1)
    {
        ITstream rkCoeffInput = multiStage.lookup("RKCoeff");

        forAll(beta,betai)
        {
            rkCoeffInput >> beta[betai];
        }
    }
}

surfaceScalarField meshPhi
(
    IOobject
    (
        "meshPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("meshPhi", dimensionSet(0,3,-1,0,0,0,0), 0.0)
);
#include "createFvOptions.H"
