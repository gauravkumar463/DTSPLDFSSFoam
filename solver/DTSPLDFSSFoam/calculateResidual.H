tresp = gMax((mag(rhoFlux)/trSubDeltaT.ref())().primitiveField());
const Vector<label>& directions = mesh.geometricD();
for (direction dir=0; dir<directions.nComponents; dir++)
{
    if (directions[dir] != -1)
    {
    	tresU.component(dir) = gMax((mag(rhoUFlux.component(dir))/trSubDeltaT.ref())().primitiveField());
	}
}
trese = gMax((mag(rhoEFlux)/trSubDeltaT.ref())().primitiveField());

tresp1 = Foam::sqrt(gAverage(sqr((mag(rhoFlux)/trSubDeltaT.ref())().primitiveField())));
for (direction dir=0; dir<directions.nComponents; dir++)
{
    if (directions[dir] != -1)
    {
    	tresU1.component(dir) = Foam::sqrt(gAverage(sqr((mag(rhoUFlux.component(dir))/trSubDeltaT.ref())().primitiveField())));
	}
}
trese1 = Foam::sqrt(gAverage(sqr((mag(rhoEFlux)/trSubDeltaT.ref())().primitiveField())));
