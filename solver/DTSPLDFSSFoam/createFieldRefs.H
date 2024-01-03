const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
//const volScalarField gamma = thermo.gamma();
//const volScalarField Cp = thermo.Cp();
//volScalarField e = thermo.he();


dimensionedScalar Pinf("Pinf", dimPressure, 0.0);
residualConvergence.lookup(Pinf.name()) >> Pinf.value();
dimensionedScalar Tinf("Tinf", dimTemperature, 0.0);
residualConvergence.lookup(Tinf.name()) >> Tinf.value();

volScalarField dp
(
    IOobject
    (
        "P_gauge",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p-Pinf
);

volScalarField dT
(
    IOobject
    (
        "T_gauge",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    T-Tinf
);


LDFSSFlux LDFSS;

volScalarField dh
(
    IOobject
    (
        "dh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
mesh,
dimensionedScalar("dh", dimLength, 0.0)
);

dh.ref() = pow(mesh.V(),1.0/3.0);

if(mesh.nGeometricD() == 2){
//if(mesh.nSolutionD() == 2){
const Vector<label>& directions = mesh.geometricD();

scalar thickness = 0.0;
for (direction dir=0; dir<directions.nComponents; dir++)
{
    if (directions[dir] == -1)
    {
        thickness = mesh.bounds().span()[dir];
        break;
    }
}
dh.primitiveFieldRef() = sqrt(mesh.V()/thickness);
}

volScalarField tp = p;
volVectorField tU = U;
volScalarField te = e;

volVectorField ddtRhoU
(
    IOobject
    (
        "ddtRhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("ddtRhoU", dimensionSet(1,-2,-2,0,0,0,0), vector::zero),
    "zeroGradient"
);
