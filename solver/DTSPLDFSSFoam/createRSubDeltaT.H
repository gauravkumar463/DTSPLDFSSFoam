tmp<volScalarField> trSubDeltaT;

if (LTS)
{

    trSubDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rSubDeltaTName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
}
