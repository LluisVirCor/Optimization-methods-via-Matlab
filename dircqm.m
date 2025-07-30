% direccions conjugades per a funcions quadratiques   
%                        requereix 'dds' amb "Q","c" i "eps"
% Inicialitzacions
kse=fopen('dircqm.res','w');
n=length(c);
x0=zeros(n,1);
x=x0;
g=Q*x+c;
ngrq=g'*g
%                            Iteracions ("n" com a molt)
ite=0;
while ((ite<n) & (ngrq>eps))
   ite=ite+1
%                            gener. de direccio Q-conj. per met. Gramm-Schmidt
   p=zeros(n,1);
   p(ite)=1;
   d=p;
   if (ite>1)
      for i=1:(ite-1)
         di=DGS(:,i);
         qdi=Q*di;
         d=d-((p'*qdi)/(di'*qdi))*di;
      end
   end
   DGS(:,ite)=d;
%                            iteracio
   qd=Q*d;
   alf=-((g'*d)/(d'*qd));
   x=x+alf*d;
   g=g+alf*qd;
   ngrq=g'*g
%                             escriptura a arxiu
   fprintf(kse,' ite= %i ngrq= %g\n',ite,ngrq);
   fprintf(kse,' x:');
   for i=1:n
      fprintf(kse,' %g',x(i)); 
   end
   fprintf(kse,'\n');
end
fclose(kse);
