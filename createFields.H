Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading Carbon electrode porosity\n" << endl;

dimensionedScalar epsilon("epsilon", dimless, transportProperties);
dimensionedScalar epsilon_0("epsilon_0", dimless, transportProperties);
dimensionedScalar alpha_2a("alpha_2a",dimless,transportProperties);
dimensionedScalar alpha_2c("alpha_2c",dimless,transportProperties);
dimensionedScalar k2("k2",dimensionSet(0,1,-1,0,0,0,0),transportProperties);
dimensionedScalar Ae("Ae",dimensionSet(0,-1,0,0,0,0,0),transportProperties);
dimensionedScalar a("a",pow(epsilon/epsilon_0,0.75)*Ae);
dimensionedScalar D_II("D_II",dimensionSet(0,2,-1,0,0,0,0),transportProperties);
dimensionedScalar D_IIeff("D_IIeff",D_II*pow(epsilon,1.5));
dimensionedScalar D_III("D_III",dimensionSet(0,2,-1,0,0,0,0),transportProperties);
dimensionedScalar D_IIIeff("D_IIIeff",D_III*pow(epsilon,1.5));

dimensionedScalar T("T",dimTemperature,transportProperties);
dimensionedScalar U0_2("U0_2",dimensionSet(1,2,-3,0,0,-1,0),transportProperties); // Volts

dimensionedScalar sigma_s("sigma_s",dimensionSet(-1,-3,3,0,0,2,0),transportProperties); // Siemens/m
dimensionedScalar sigma_seff("sigma_seff",pow(1.0-epsilon,1.5)*sigma_s);
// We should check weather the following variables should be dimensionaedScalar or volScalarField
dimensionedScalar F("F",dimensionSet(0,0, 1, 0,-1,1,0),96485.3321233100184);
dimensionedScalar R("R",dimensionSet(1,2,-2,-1,-1,0,0),8.31446261815324);

// i0_2=F*k2*pow(V_II,alpha_2c)*pow(V_III,alpha_2a);




Info<< "Reading field V_II\n" << endl;

volScalarField V_II
(
    IOobject
    (
        "V_II",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field V_III\n" << endl;

volScalarField V_III
(
    IOobject
    (
        "V_III",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField phi_S
(
    IOobject
    (
        "phi_S",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField phi_L
(
    IOobject
    (
        "phi_L",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating field U_2" << endl;
volScalarField U_2
(
    IOobject
    (
        "U_2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U0_2+(R*T/F)*log(V_III/V_II)
);

Info << "Creating field eta_2" << endl;
volScalarField eta_2
(
    IOobject
    (
        "eta_2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi_S-phi_L-U_2
);

volScalarField i0_2
(
    IOobject
    (
        "i0_2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    // mesh,
    F*k2*pow(V_II,alpha_2c)*pow(V_III,alpha_2a)
);

volScalarField j_2
(
    IOobject
    (
        "j_2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    i0_2*(exp(-alpha_2c*F*eta_2/(R*T))-exp(alpha_2a*F*eta_2/(R*T)))
);

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


Info << "Creating field k_Leff" << endl;
volScalarField k_Leff
(
    IOobject
    (
        "k_Leff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pow(F,2)/(R*T)*(D_IIeff*V_II+D_IIIeff*V_III) // Ojo que falta z_i**2
);
#include "createPhi.H"

// #include "createFvOptions.H"
