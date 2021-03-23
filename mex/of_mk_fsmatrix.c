
/* ------------------------------------------------------------------------*/
/* Step0 declarations at the beginning of 'MatrixIntegration' */

/* ------------------------------------------------------------------------*/
/* STEP1 prepare for strategy */ 
#if  (MatrixIntegrationStep==1)
/* Follower pressure (fluid/structure coupling) */
else if (!strcmp("fs_matrix",CAM)) { 
       StrategyType=4;
     if (pointg[4]==5) {
       Ndof = 4*Nnode; GF.NBe=Ndof; GF.NdefE=Ndof*Ndef;
     }
}
#endif

/* ------------------------------------------------------------------------*/
/* STEP2 switch for strategy implementation  */ 
#if  (MatrixIntegrationStep==2)
case 4: {

  /* compute matrix and rhs for follower pressure, see of_mk_pre.c */
  F_pressure(GF,ECp,Nnode,Nw,point);

}
break; 

#endif

/* ------------------------------------------------------------------------*/
/* STEP3 possible cleanup if needed */

