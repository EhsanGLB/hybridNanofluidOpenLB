/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.1                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
nanoFluid4OLB
{
    version     1.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//name rho Cp kappa mu beta;
water 997.1 4179 0.613 0.0009087 0.00021;
SiO2 2220 745 1.4 0.0 0.0000055;
MWCNT 1600 796 3000 0.008 0.0000085;


lambdanp1	3.0;
lambdanp2	6.0;
phinp1		0.02;
phinp2		0.02;


Tcold		283.15;
Thot		293.15;
g		9.81;
Ra		1000000.0;


pointNum	200;
endTime		200;
residual	0.01;
