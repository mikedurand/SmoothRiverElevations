

reAttach=true;
p.distmax=500; 
p.ncut=20;
p.StdMax=2.5;
p.pct=25;
p.dx=.1;
p.x=65:p.dx:90;
p.N=length(p.x);
p.HrefCut=1;

Href{1}=nan(1,p.N); Href{2}=nan(1,p.N);
[Est]=ProcessData(c,s,reAttach,p,Href,'DataName','LP');

%%
reAttach=false;
Href{1}=Est{1}.Hc; Href{2}=Est{2}.Hc;
[Est,Data]=ProcessData(c,s,reAttach,p,Href,'DataName','SLM');

