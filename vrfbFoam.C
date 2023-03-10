/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vrfbFoam

Description
    Vanadium redox flow battery solver.

    This solver aims at solve a half battery

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Vanadium Redox Flow Battery solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting the VRFB solver\n" << endl;

    // #include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        i0_2=F*k2*pow(V_II,alpha_2c)*pow(V_III,alpha_2a);
        k_Leff=pow(F,2.0)/(R*T)*(D_IIeff*V_II+D_IIIeff*V_III); // Ojo que falta z_i**2
        U_2=U0_2+(R*T/F)*log(V_III/V_II);
        U_2.relax();
        eta_2=phi_S-phi_L-U_2;
        // Info << alpha_2a << endl;
        // Info << F << endl;
        // Info << R << endl;
        // Info << T << endl;
        // Info << "eta_2" << endl;
        // Info << eta_2 << endl;
        // Info << "U_2" << endl;
        Info << (eta_2*F)/(R*T) << endl;
        //j_2=i0_2*(exp(-alpha_2c*F*eta_2/(R*T)));//-
        j_2=i0_2*(exp(alpha_2a*F*eta_2/(R*T)));
        Info << j_2 << endl;
        // j_2=i0_2*(exp(alpha_2a*F*eta_2/(R*T)));


        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix V_IIEqn
            (
                fvm::div(phi, V_II)
              - fvm::laplacian(D_IIeff, V_II)
             ==
                a*j_2/F
            );

            V_IIEqn.relax();
            // fvOptions.constrain(TEqn);
            V_IIEqn.solve();
            // fvOptions.correct(T);

            fvScalarMatrix V_IIIEqn
            (
                fvm::div(phi, V_III)
              - fvm::laplacian(D_IIIeff, V_III)
             ==
                -a*j_2/F
            );

            V_IIIEqn.relax();
            // fvOptions.constrain(TEqn);
            V_IIIEqn.solve();
            // fvOptions.correct(T);

            fvScalarMatrix phi_S_Eqn
            (
              - fvm::laplacian(sigma_seff, phi_S)
             ==
                a*j_2
            );

            phi_S_Eqn.relax();
            // fvOptions.constrain(TEqn);
            phi_S_Eqn.solve();
            // fvOptions.correct(T);

            fvScalarMatrix phi_L_Eqn
            (
              - fvm::laplacian(k_Leff, phi_L)
             ==
                -a*j_2
            );

            phi_L_Eqn.relax();
            // fvOptions.constrain(TEqn);
            phi_L_Eqn.solve();
            // fvOptions.correct(T);

        }

        // The following is to delete negative values of V_II and V_III
        V_II.max(SMALL);
        V_III.max(SMALL);

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
