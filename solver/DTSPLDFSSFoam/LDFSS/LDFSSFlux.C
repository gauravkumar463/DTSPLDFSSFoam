/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    AUSMplusPreFlux

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "LDFSSFlux.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //+


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LDFSSFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    scalar& amaxSf,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
    const scalar TLeft,
    const scalar TRight,
    const scalar Ur2Sf,
    const scalar gamma,
    const scalar Cp,
    const vector Sf,
    const scalar phi,
    const scalar Pinf,
    const scalar Tinf
) const
{
    // Info << rhoFlux << " " << rhoUFlux << " " << rhoEFlux << " "
    //      << amaxSf << " " << pLeft << " " << pRight << " "
    //      << ULeft << " " << URight << " " << TLeft << " "
    //      << TRight << " " << Ur2Sf << " " << gamma << " "
    //      << Cp << " " << Sf << endl;
    scalar magSf = mag(Sf);
    vector n = Sf/magSf;

    scalar R = Cp*(gamma-1.0)/gamma;

    scalar rhoLeft = (Pinf+pLeft)/(R*(Tinf+TLeft));
    scalar rhoRight = (Pinf+pRight)/(R*(Tinf+TRight));
    scalar rhoHalf = 0.5*(rhoLeft+rhoRight);

    scalar aLeft = sqrt((gamma-1)*Cp*(Tinf+TLeft));
    scalar aRight = sqrt((gamma-1)*Cp*(Tinf+TRight));
    scalar aHalf = sqrt(0.5*(sqr(aLeft)+sqr(aRight)));
    //scalar aHalf = 0.5*(aLeft+aRight);
    //scalar aHalf = sqrt(aLeft*aRight);

    scalar qLeft = (ULeft & n)-(phi/magSf);
    scalar qRight = (URight & n)-(phi/magSf);

    scalar MLeft = qLeft/aHalf;
    scalar MRight = qRight/aHalf;
    scalar Mavg = sqrt(0.5*(sqr(MLeft)+sqr(MRight)));
    //scalar Mavg = 0.5*(MLeft+MRight);
    //scalar Mavg = sqrt(MLeft*MRight);
    scalar Mr2Half = Ur2Sf/sqr(aHalf);

    scalar fHalf =  sqrt((sqr(1-Mr2Half)*sqr(Mavg))+(4*Mr2Half))/(1+Mr2Half);
    scalar aTildeHalf = fHalf*aHalf;

    MLeft = MLeft/fHalf;
    MRight = MRight/fHalf;

    amaxSf = 0.5*aTildeHalf*(1+Mr2Half)*(Mavg+1)*magSf;
    //amaxSf = (0.5*(1.0+Mr2Half))*((0.5*(qLeft+qRight)) + aTildeHalf)*magSf;

    scalar alphaLeft = 0.5*(1.0 + sign(MLeft));
    scalar alphaRight = 0.5*(1.0 - sign(MRight));

    scalar betaLeft = -max(0.0,1-(int)mag(MLeft));
    scalar betaRight = -max(0.0,1-(int)mag(MRight));

    //scalar M1Left = 0.5*(MLeft+mag(MLeft));
    //scalar M1Right = 0.5*(MRight-mag(MRight));

    scalar M4Left = 0.25*sqr(MLeft+1) + (1.0/8.0)*sqr(sqr(MLeft)-1);
    scalar M4Right= -0.25*sqr(MRight-1) - (1.0/8.0)*sqr(sqr(MRight)-1);

    scalar Mplus = (alphaLeft*(1+betaLeft)*MLeft) - (betaLeft*M4Left);
    scalar Mminus = (alphaRight*(1+betaRight)*MRight) - (betaRight*M4Right);

    //scalar MHalf = 0.5*(Mplus-Mminus - (alphaLeft*MLeft) + (alphaRight*MRight));
    scalar MHalf = 0.25*betaLeft*betaRight*sqr(sqrt(0.5*(sqr(MLeft)+sqr(MRight)))-1.0);
    //scalar MHalf = 0.25*betaLeft*betaRight*sqr(0.5*(MLeft+MRight)-1.0);

    scalar UPlus = aTildeHalf*(Mplus - (MHalf*(1.0-((pLeft-pRight)/(2*rhoLeft*Ur2Sf)))));
    scalar Uminus = aTildeHalf*(Mminus + (MHalf*(1.0+((pLeft-pRight)/(2*rhoRight*Ur2Sf)))));

    scalar M5Left = 0.5*(1.0+MLeft);
    scalar M5Right = 0.5*(1.0-MRight);
    //scalar M5Left = 0.25*sqr(MLeft+1)*(2-MLeft);// + (3.0/16.0)*MLeft*sqr(sqr(MLeft)-1);
    //scalar M5Right = -0.25*sqr(MRight-1)*(2+MRight);// - (3.0/16.0)*MRight*sqr(sqr(MRight)-1);

    scalar p1Left = alphaLeft*(1+betaLeft) - betaLeft*M5Left;
    scalar p1Right = alphaRight*(1+betaRight) - betaRight*M5Right;

    scalar pHalf = Pinf + 0.5*(pLeft+pRight) + 0.5*(p1Left-p1Right)*(pLeft-pRight) 
                                + rhoHalf*Ur2Sf*(p1Left+p1Right-1.0);

    scalar HLeft = (Cp/gamma)*(Tinf+TLeft) + 0.5*magSqr(ULeft) + pHalf/rhoLeft;
    scalar HRight = (Cp/gamma)*(Tinf+TRight) + 0.5*magSqr(URight) + pHalf/rhoRight;

    rhoFlux = magSf*(rhoLeft*UPlus + rhoRight*Uminus);

    rhoUFlux = magSf*(rhoLeft*ULeft*UPlus + rhoRight*URight*Uminus) + pHalf*Sf;

    rhoEFlux = magSf*(rhoLeft*HLeft*UPlus + rhoRight*HRight*Uminus);
}

// ************************************************************************* //
