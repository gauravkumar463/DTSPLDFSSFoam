// calculate square of speed of sound
volScalarField c2 = (thermo.gamma()-1)*thermo.Cp()*T;

volScalarField U2 = magSqr(U);

dimensionedScalar Uinf("Uinf", dimVelocity, 0.0);
residualConvergence.lookup(Uinf.name()) >> Uinf.value();

if (mesh.moving())
{
    dh.ref() = pow(mesh.V(),1.0/3.0);

    if(mesh.nGeometricD() == 2){
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

    meshPhi = mesh.phi();
}

if(preconditioning == false){
    Ur2.primitiveFieldRef() = c2.primitiveField();
}
else{
    if(!inviscid == true){
    Ur2.primitiveFieldRef() = min
                      (
                      c2.primitiveField(), max
                              (sqr(Uinf.value()), max
                                     (U2.primitiveField(), max
                                          (1e-6*mag(fvc::grad(dp))().primitiveField()*dh.primitiveField()/rho.primitiveField(), sqr(max
                                              (
                                                turbulence->nuEff()/(dh.primitiveField()),
                                                turbulence->alphaEff()/(rho*dh.primitiveField())
                                              ))
                                          )
                                    )
                              )
                      );
    }
    else{
        Ur2.primitiveFieldRef() = min
                      (
                      c2.primitiveField(), max
                              (sqr(Uinf.value()), max
                                     (U2.primitiveField(), 1e-6*mag(fvc::grad(dp))().primitiveField()*dh.primitiveField()/rho.primitiveField())
                              )
                      );
    }
}
Ur2.correctBoundaryConditions();
// --- Directed interpolation of primitive fields onto faces

surfaceScalarField T_pos(interpolate(dT, pos, T.name()));
surfaceScalarField T_neg(interpolate(dT, neg, T.name()));

surfaceScalarField p_pos(interpolate(dp, pos, p.name()));
surfaceScalarField p_neg(interpolate(dp, neg, p.name()));

surfaceVectorField U_pos(interpolate(U, pos, U.name()));
surfaceVectorField U_neg(interpolate(U, neg, U.name()));

surfaceScalarField Ur2Sf(sqr(linearInterpolate(sqrt(Ur2))));

surfaceScalarField gammaSf(linearInterpolate(thermo.gamma()));
surfaceScalarField CpSf(linearInterpolate(thermo.Cp()));

forAll(phi,faceI){
    LDFSS.evaluateFlux
    (
        phi[faceI],
        phiUp[faceI],
        phiEp[faceI],
        amaxSf[faceI],
        p_pos[faceI],
        p_neg[faceI],
        U_pos[faceI],
        U_neg[faceI],
        T_pos[faceI],
        T_neg[faceI],
        Ur2Sf[faceI],
        gammaSf[faceI],
        CpSf[faceI],
        mesh.Sf()[faceI],
        meshPhi[faceI],
        Pinf.value(),
        Tinf.value()
    );
}

surfaceScalarField::Boundary& bfield = phi.boundaryFieldRef();

forAll(bfield, patchi)
{
    const scalarField& currPatch = phi.boundaryField()[patchi];

    fvsPatchScalarField& phip = phi.boundaryFieldRef()[patchi];
    fvsPatchVectorField& phiUpp = phiUp.boundaryFieldRef()[patchi];
    fvsPatchScalarField& phiEpp = phiEp.boundaryFieldRef()[patchi];
    fvsPatchScalarField& amaxSfp = amaxSf.boundaryFieldRef()[patchi];
    const fvPatchScalarField& pp = dp.boundaryField()[patchi];
    const fvPatchVectorField& Up = U.boundaryField()[patchi];
    const fvPatchScalarField& Tp = dT.boundaryField()[patchi];
    const fvPatchScalarField& Ur2p = Ur2.boundaryField()[patchi];
    const fvsPatchScalarField& gammap = gammaSf.boundaryField()[patchi];
    const fvsPatchScalarField& Cpp = CpSf.boundaryField()[patchi];
    const fvsPatchVectorField& Sfp = mesh.Sf().boundaryField()[patchi];
    const fvsPatchScalarField& meshPhip = meshPhi.boundaryField()[patchi];
    
    if (bfield[patchi].coupled())
    {
        const scalarField pcp = dp.boundaryField()[patchi].patchInternalField();
        const vectorField Ucp = U.boundaryField()[patchi].patchInternalField();
        const scalarField Tcp = dT.boundaryField()[patchi].patchInternalField();
        const scalarField pcnb = dp.boundaryField()[patchi].patchNeighbourField();
        const vectorField Ucnb = U.boundaryField()[patchi].patchNeighbourField();
        const scalarField Tcnb = dT.boundaryField()[patchi].patchNeighbourField();

        forAll(currPatch,faceI)
        {
            LDFSS.evaluateFlux
            (
                phip[faceI],
                phiUpp[faceI],
                phiEpp[faceI],
                amaxSfp[faceI],
                pcp[faceI],
                pcnb[faceI],
                Ucp[faceI],
                Ucnb[faceI],
                Tcp[faceI],
                Tcnb[faceI],
                Ur2p[faceI],
                gammap[faceI],
                Cpp[faceI],
                Sfp[faceI],
                meshPhip[faceI],
                Pinf.value(),
                Tinf.value()
            );
        }
    }
    else
    {
        forAll(currPatch,faceI)
        {
            LDFSS.evaluateFlux
            (
                phip[faceI],
                phiUpp[faceI],
                phiEpp[faceI],
                amaxSfp[faceI],
                pp[faceI],
                pp[faceI],
                Up[faceI],
                Up[faceI],
                Tp[faceI],
                Tp[faceI],
                Ur2p[faceI],
                gammap[faceI],
                Cpp[faceI],
                Sfp[faceI],
                meshPhip[faceI],
                Pinf.value(),
                Tinf.value()
            );
        }

    }
}
