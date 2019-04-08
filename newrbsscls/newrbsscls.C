/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    int count=0;
    if (!LTS)
    {
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    }

    #include "mappingPsi.H" //scls
    #include "solveLSFunction.H"  //scls
    #include "calcNewCurvature.H" //scls

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //========================================================================
    // initialse the file that will dump the total mass
    
    autoPtr<OFstream> massFilePtr;
    
    scalar totalMass(0.0);
    scalar totalMass0(0.0);


   const scalarField& V = mesh.V() ;
    totalMass0 = gSum(rho*V);
    
    
    //========================================================================



    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        

     if (LTS)
        {
      #include "setRDeltaT.H"
         }
        else
        {

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
          }
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

	H.storeOldTime(); //scls

        // --- Pressure-velocity PIMPLE corrector loop
        //while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"



        #include "mappingPsi.H" //scls
        #include "solveLSFunction.H" //scls
        #include "calcNewCurvature.H" //scls
        #include "updateFlux.H" //scls

            mixture.correct();

    forAll(mesh.cells(),celli)
    {
        if (C[celli] > scalar(1e3)) 
          {
             count =  1;
          }
        
     }


	if (count==0)
	  {
        #include "UEqn.H"

        #include "pEqn.H"
        }
	if (count==1)
	  {
        #include "U1Eqn.H"

        #include "p1Eqn.H"
        }


        count=0;
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

	// reInitialise the alpha equation 
	if (runTime.outputTime())
	  {
	    Info<<"Overwriting alpha" << nl << endl;
	    alpha1 = H;
            volScalarField& alpha10 = const_cast<volScalarField&>(alpha1.oldTime());
	    alpha10 = H.oldTime();
	  }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
