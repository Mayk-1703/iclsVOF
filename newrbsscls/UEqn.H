Info<< "ueq " <<  nl << endl;     
MRF.correctBoundaryVelocity(U);
//Info<< "ueq " <<  nl << endl;         
fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvc::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
   //  ==
      //  fvOptions(rho, U)
    );
//Info<< "ueq " <<  nl << endl;     
    UEqn.relax();
//Info<< "ueq " <<  nl << endl;     
    //fvOptions.constrain(UEqn);

    //if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (

                    sigma*fvc::snGrad(psi)*fvc::interpolate(C)*fvc::interpolate(delta)
		  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        //fvOptions.correct(U);
    }
//Info<< "ueqend " <<  nl << endl;  
