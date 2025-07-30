% gradient conjugat per a funcions quadratiques   
%                        requereix 'dds' amb "Q","c" i "eps"
% Inicialitzacions
kse=fopen('gracqm.res','w');
n=length(c);
x0=zeros(n,1);
x=x0;
g=Q*x+c;
ngrq=g'*g
d=-g;
%                            Iteracions ("n" com a molt)
ite=0;
while ((ite<n) & (ngrq>eps))
   ite=ite+1
   qd=Q*d;
   alf=ngrq/(d'*qd);
   x=x+alf*d;
   g=g+alf*qd;
   ngacq=g'*g                % norma del gradient al quadrat actual
    
   DGC(:,ite)=d;
   beta=(g'*qd)/(d'*qd);
   d=-g+beta*d;
   ngrq=ngacq;               % norma del gradient al quadrat anterior
%                             escriptura a arxiu
    fprintf(kse,' ite= %i ngrq= %g\n',ite,ngrq);
    fprintf(kse,' x:');
    for i=1:n
       fprintf(kse,' %g',x(i)); 
    end
    fprintf(kse,'\n');
end
% si hi ha hagut menys de "n" iteracions, mirem els valors propis de "Q"
if (ite<n)
   vaps=eig(Q)
end
fclose(kse);
