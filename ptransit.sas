
/****************************************************************************;
Macro name: %ptransit

Function: Returns predicted transition probabilities for a proportional hazards model
          for recurrent events by using the methodology proposed by Andersen, Hansen 
          and Keiding (1991). Given a data set in the counting process style of input
          with multiple observations for each individual, the macro fits the model and 
          provides estimates and graphics for subjects from a new data set.  

		  The models considered are the independent increments (AG), 
          conditional (PWP) and marginal (WLW) proposed respectively by
          Andersen and Gill (1982), Prentice, Williams and Peterson (1981) and 
          Wei, Lin and Weissfeld (1989).	

		  Some events are immediately reversible and the interest is only in the
          probability of recurrence. Others are not immediatelly reversible, and 
          the probability of recovery can be of interest. The macro deals with both 
          and estimates one or two types of transition probabilities depending on 
          the choice in the option "RECOV".
          

Developer:    A. T. Paes
              Laboratory of Epidemiology and Statistics
              Institute Dante Pazzaneze of Cardiology   
              Sao Paulo, Brazil

Date:  January 25, 2002

Call statement:

           %ptransit (data=,time1=,time2,event=,xvars=,
                      id=,model=,new=,recov=,strata=,transit=)

           Parameters required:


              data = is the data set used to fit the model. The style of data
                     input must follow the counting process formulation, where each
                     subject is represented by multiple observations. In this style, as the
                     data for an individual are displayed in several rows, the values of
                     time-independent covariates must be repeated in all of them.

              time1 and time2 = are the variables representing the time interval during which each
                                observation is at risk . TIME1 is the left endpoint and TIME2 is the right endpoint of the risk interval. For
                                studies with only one type of transition, even when there is no discontinuity of the risk intervals, the data
                                set must include a variable corresponding to TIME1 where all values are 0.

              event = is the variable that indicates whether the
                      transition of interest occurred at the end of the risk interval. It
                      assumes the values 0 and 1.

              xvars = is the list of covariates to be used in the regression model.

              id = is the subject identification. Since the
                   counting process style of input is used, the ID is repeated as many
                   times as the observations for each individual.

              model= defines which model must be used to fit the
                     parameters: AG, PWP or WLW. By default, the AG model is assumed.

              new = is the data set containing the covariate values
                    for which the transition probabilities are to be estimated. It may
                    contain hypothetical observations, covariate values from the original
                    data set or any meaningful values one is interested on. This data set
                    encloses the values of covariates and each row corresponds to each
                    individual.  The order of the covariates columns must be the same as
                    the original data set.

              recov = indicates whether the ``diseased $\rightarrow$ healthy'' transition is present, that
                      is, defines the number of transitions of interest. The option RECOV=Y means that there are two types of
                      transitions and this instructs the program to create a new data set with {\it type-specific} covariates which
                      assume the original value for the specified transition and 0 otherwise. By default, it is assumed to be N (one
                      transition).

              strata = represents the variable corresponding to the
                       event order or transition occurrences (first time, second time, third
                       time, etc.). If specified, the model assumes different baselines
                       intensities for each stratum. It is not used for the AG model.

              transit = is the name of the variable that distinguish
                        the type of transition. For example, if there are two types of
                        transition (0 -> 1 and 1 -> 0) this variable
                        must assume the value 1 for 0 -> 1 and the value 2 for 1 -> 0
                        transition. It does not exist when there is only one
                        type of transition.

Output: The macro prints the PHREG output used to fit the Cox model and a
        table including the predicted transition probabilities for each 
        individual from the new data set at all the observed transition times. It also
        includes a graphic with these estimates in which each line represents 
        an individual. In case of two types of transition, two graphics are provided.


References:


         Andersen, P.K. and Hansen, L.S. and Keiding, N. (1991).
         "Non-and Semi-parametric Estimation of Transition
         Probabilities from Censored Observation of a Non-homogeneous {M}arkov Process",
         Scandinavian Journal of Statistics,18,153-167.

         Andersen, P.K. and Gill, R.D. (1982)."Cox's Regression Model 
         for Counting Processes: A Large Sample Study", 
         The Annals of Statistics, 10, 1100-1120.

         Prentice, R.L. and Williams, B.J. and Peterson, A.V. (1981).
         "On the Regression Analysis of Multivariate Failure Time Data",
         Biometrika, 68,373-379.
         
         Wei, L.J. and Lin, D.Y. and Weissfeld, L. (1989).
         "Regression Analysis of Multivariate Incomplete Failure Time Data 
         by Modeling Marginal Distributions",JASA, 84, 1065-1073.
         
         Therneau, T.M. and Grambsch, P.M. (2000).
         "Moddeling Survival Data - Extending the Cox Model",
         Springer-Verlag", New York.
      
***************************************************************************;
     */


%macro ptransit(data=,time1=,time2=,event=,xvars=,
id=,model=,new=,recov=,strata=,transit=);


/* count the number of covariates */

%let nxs=0;
%do i=1 %to 50;
   %if %scan(&xvars,&i)= %then %goto done;
   %let nxs=%eval(&nxs+1);
%end;
%done: %put &nxs;


/* create the name of specific covariates (only for two transitions) */
/* Create a new data set including the specific covariates */
/* identify distinct subjects */

%if %upcase(&recov)=Y %then %do;
   %local i;
   %let names= ;
   %do j=1 %to 2;
    %do i=1 %to &nxs;
     %local var names;
     %let var=%scan(&xvars,&i)&j;
     %let names=&names &var;
    %end;
   %end;
   %put &names;

   data data1;
   set &data;
   by &id;
   %local k;
   %local var var2;
   %let k=1;
   %do j=1 %to 2;
     %do i=1 %to &nxs;
       %let var=%scan(&names,&k);
       %let var2=%scan(&xvars,&i);
       &var = &var2 * (&transit=&j);
       %let k=%eval(&k+1);
     %end;
   %end;
   primind=first.&id;
   run;
%end;

%else %do;
  data data1;
  set &data;
  by &id;
  primind=first.&id;
  run;
%end;


/* Fit the Cox model */


%if (%upcase(&model)=WLW) %then %do;
 proc phreg data=data1 out=estim COVSANDWICH(AGGREGATE);
 %if %upcase(&recov)=Y %then %do;
  model (&time1,&time2)*&event(0)= &names
 %end;
 %else %do;
  model (&time1,&time2)*&event(0)= &xvars
 %end;
 /risklimits alpha=0.05;
 id &id;
 strata &strata &transit;
 run;
%end;


%else %do; 
  proc phreg data=data1 out=estim;
  %if %upcase(&recov)=Y %then %do;
   model (&time1,&time2)*&event(0)= &names
  %end;
  %else %do;
   model (&time1,&time2)*&event(0)= &xvars
  %end;
  /risklimits alpha=0.05;
  id &id;
    strata &strata &transit;
  run;
%end; 



/* Sort the data set by the observed transition times */

proc sort data=data1 out=data2;
by &time2;
run;


/* Create the variable with the first observation at each transition time and
variables that indicates the ocurrences of the transitions at each observation */


data tempos;
set data2;
%if %upcase(&recov)=Y %then %do;
keep &id &transit &event &time2 primt N01 N10;
%end;
%else %do;
keep &id &event &time2 primt N01;
%end;
by &time2;
primt=first.&time2;
%if %upcase(&recov)=Y %then %do;
N01=1*(&transit=1)*(&event=1);
N10=1*(&transit=2)*(&event=1);
%end;
%else %do;
N01=1*(&event=1);
%end;
run;



/****************************************************************/
/*                                                              */
/* Compute the vector TIME with the transition times            */
/* DN01 (number of 0-1 transitions at each time                 */ 
/****************************************************************/



proc iml;
use tempos;

read all var{&time2} where (primt=1) into tempo; /* vector of observed transition times     */
ntemp=nrow(tempo);                               /* number of distinct observed times */
read all var{&time2} into temp;                  /* vector with all the times (duplicate) */
ntemps=nrow(temp);                               /* total number of observations */
read all var{N01} into N01;       /* vector with variables that indicate the 0-1 transition */
%if %upcase(&recov)=Y %then %do;
read all var{N10} into N10;  /*    "                 "                       1-0 transition */
%end;


DN01=J(ntemp,1,0);
sum01=0;
%if %upcase(&recov)=Y %then %do;
DN10=J(ntemp,1,0);
sum10=0;
%end;




do T1=1 to ntemp;
    do T2=1 to ntemps;
    if (temp[T2]=tempo[T1]) then do;
	   sum01=sum01+N01[T2];
       %if %upcase(&recov)=Y %then %do;
       sum10=sum10+N10[T2];
	   %end;
    end;
   end;

   DN01[T1]=sum01;   /* vector with the number of transitions at each time */
   sum01=0;
   %if %upcase(&recov)=Y %then %do;
   DN10[T1]=sum10;
   sum10=0;
   %end;
end;


/*************************************************************************************/
/*                                                                                   */
/* Create the matrices Y01 e Y10 that indicate at each time whether the individuals   */
/* are at risk for the 0-1 ou 1-0 transitions.                                      */
/* Obs: dimension of the matrices: (number of individuals) x (number of distinct times) */
/*                                                                                   */
/*************************************************************************************/


use data1;
read all var{&id} where (primind=1) into id;
read all var{&id} into regist;
read all var{&time2} into tempo2;
%if %upcase(&recov)=Y %then %do;
read all var{&transit} into tipotran;
%end;

nid=nrow(id);
nlinh=nrow(regist);
Y01=J(nid,ntemp,1);
%if %upcase(&recov)=Y %then %do;
Y10=J(nid,ntemp,0);
%end;


do count=1 to nlinh;
   I1=loc(id=regist[count]);
   K1=loc(tempo=tempo2[count]);
   if count < nlinh then do;
         I2=loc(id=regist[count+1]);
		 %if %upcase(&recov)=Y %then %do;
         K2=loc(tempo=tempo2[count+1]);
		 %end;
   end;

   %if %upcase(&recov)=Y %then %do;
   if I2 = I1 then do;
       if (tipotran[count]=1) then do;
          do M=K1+1 to K2;
            Y01[I1,M]=0;
            Y10[I1,M]=1;
            end;
       end;
   end;
   else do;
     do M=K1+1 to ntemp;
     Y01[I1,M]=0;
     Y10[I1,M]=0;
     end;
   end;
%end;
%else %do;
  if I2 ^= I1 then do;
      do M=K1+1 to ntemp;
      Y01[I1,M]=0;
      end;
  end;
%end;

end;

sumY01=J(ntemp,1,0);
J0=J(ntemp,1,1);
soma1=0;
%if %upcase(&recov)=Y %then %do;
sumY10=J(ntemp,1,0);
J1=J(ntemp,1,0);
soma2=0;
%end;

do aux1=1 to ntemp;
   do aux2=1 to nid;
      soma1=soma1+Y01[aux2,aux1];
	  %if %upcase(&recov)=Y %then %do;
      soma2=soma2+Y10[aux2,aux1];
	  %end;
   end;
   sumY01[aux1]=soma1;
   %if %upcase(&recov)=Y %then %do;
   sumY10[aux1]=soma2;
   %end;
   if sumY01[aux1]=0 then J0[aux1]=0;
   soma1=0;
   %if %upcase(&recov)=Y %then %do;
   if sumY10[aux1]>0 then J1[aux1]=1;
   soma2=0;
   %end;
end;


/* Compute the linear predictors */

%if %upcase(&recov)=Y %then %do;
 read all var{&xvars} where (primind=1) into VAR;
 p=ncol(var);
 zeros=J(nid,p,0);
 Z1= VAR ||  zeros;
 Z2= zeros ||  VAR;
 use estim;
 read all var{&names} into b;
 ZB01=Z1*b`;
 ZB10=Z2*b`;
 S01=Y01`*exp(ZB01);
 S10=Y10`*exp(ZB10);
 do aux3=1 to ntemp;
  if S01[aux3]=0 then S01[aux3]=1;
  if S10[aux3]=0 then S10[aux3]=1;
 end;
 DA010=(J0/S01)#DN01;
 DA100=(J1/S10)#DN10; 
 use &new;
 read all var{&xvars} into Zteste;
 nid1=nrow(Zteste);
 zeros=J(nid1,&nxs,0);
 ZTESTE1=ZTESTE || zeros;
 ZTESTE2=zeros || ZTESTE;
 p010=Zteste1*b`;
 p100=Zteste2*b`;
 ep01=(exp(p010))`;
 ep10=(exp(p100))`;
 DA01= DA010 * ep01;
 DA00=-DA01;
 DA10= DA100 * ep10;
 DA11=-DA10;
 P=I(2);
 P01=J(ntemp,nid1,0);
 P10=J(ntemp,nid1,0);
 aux3=ntemp*nid1;
 PROB01=J(AUX3,1,0);
 PROB10=J(AUX3,1,0);
 AUX_T=J(AUX3,1,0);
 AUX_I=J(AUX3,1,0);
 i=1;
 j=1;
 do aux2=1 to nid1;
   do aux1=1 to ntemp;
     DA=(DA00[aux1,aux2] || DA01[aux1,aux2]) // (DA10[aux1,aux2] || DA11[aux1,aux2]);
     IDENT=I(2);
     P=P*(IDENT+DA);
     PROB01[j]=P[1,2];
	 PROB10[j]=P[2,1];
	 AUX_T[j]=TEMPO[aux1];
     AUX_I[j]=i;
	 j=j+1;
   end;
   P=I(2);
  i=i+1;
 end;
 do aux4=1 to aux3;
  if PROB01[aux4]<0 then PROB01[aux4]=-0.01;
  if PROB10[aux4]<0 then PROB10[aux4]=-0.01;
  if PROB01[aux4]>1 then PROB01[aux4]=1.01;
  if PROB10[aux4]>1 then PROB10[aux4]=1.01;
 end; 
Probs=AUX_I || AUX_T || PROB01 || PROB10 ;
pname={"Id","Time","Prob01","Prob10"};
create Probs from Probs [colname=pname];
append from Probs[colname=pname];
proc print data=Probs;
run;

proc gplot data=Probs;
plot prob01*time=id
/vaxis=0 to 1 by 0.1;
/* symbol1 interpol=stepLJ h=1 v=square  c=black;
 symbol2 interpol=stepLJ h=1 v=star  c=black;
 symbol3 interpol=stepLJ h=1 v=circle  c=black;
 symbol4 interpol=stepLJ h=1 v=diamond c=black;
*/
plot prob10*time=id
/vaxis=0 to 1 by 0.1; 
/* symbol1 interpol=stepLJ h=1 v=square  c=black;
 symbol2 interpol=stepLJ h=1 v=star  c=black;
 symbol3 interpol=stepLJ h=1 v=circle  c=black;
 symbol4 interpol=stepLJ h=1 v=diamond c=black;
*/
run;

%end;

%else %do;
 read all var{&xvars} where (primind=1)  into Z1;
 use estim;
 read all var{&xvars} into b;
 ZB01=Z1*b`;
 S01=Y01`*exp(ZB01);
 DA010=(J0/S01)#DN01;
 use &new;
 read all var{&xvars} into Z010;
 nid1=nrow(Z010);
 p010=Z010*b`;
 ep01=(exp(p010))`;
 DA01= DA010 * ep01;
 DA00=-DA01;
 PR01=J(ntemp,nid1,0);
 aux3=ntemp*nid1;
 PROB01=J(AUX3,1,0);
 AUX_T=J(AUX3,1,0);
 AUX_I=J(AUX3,1,0);
 P00=1;
 P01=0;
 i=1;
 j=1;
 do aux2=1 to nid1;
   do aux1=1 to ntemp;
     P00=P00*(1- DA01[aux1,aux2]);
     P01=1-P00;
     PROB01[j]=P01;
	 AUX_T[j]=TEMPO[aux1];
     AUX_I[j]=i;
	 j=j+1;
   end;
  i=i+1;
  P00=1;
  P01=0;
 end;
 Prob01=AUX_I || AUX_T || PROB01;
 pname={"Id","Time","Prob01"};
 create Prob01 from Prob01 [colname=pname];
 append from Prob01[colname=pname];
 proc print data=Prob01;
 run;

 proc gplot data=Prob01;
 plot prob01*time=id
 /vaxis=0 to 1 by 0.1;
/* symbol1 interpol=stepLJ h=1 v=square  c=black;
 symbol2 interpol=stepLJ h=1 v=diamond c=black;
 symbol3 interpol=stepLJ h=1 v=circle  c=black;
 symbol4 interpol=stepLJ h=1 v=star  c=black;
*/
run;

%end;

quit;


%mend ptransit;
