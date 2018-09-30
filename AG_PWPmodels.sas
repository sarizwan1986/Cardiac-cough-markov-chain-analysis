
PROC IMPORT OUT= WORK.COUGH 
            DATAFILE= "/folders/myfolders/sasuser.v94/COUGH.xls" 
            DBMS=XLS REPLACE;
     GETNAMES=YES;
RUN;

%include "/folders/myfolders/sasuser.v94/ptransit.sas";

%ptransit(data=COUGH,time1=time1,time2=time2,event=event,
xvars=age gender lv bh vd sur,id=id,new=NEW1,
recov=Y,transit=typetran);

%ptransit(data=COUGH,time1=time1,time2=time2,event=event,
xvars=age gender lv bh vd sur,id=id,model=PWP,new=NEW1,
recov=Y,strata=regist,transit=typetran);
