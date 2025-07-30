% algorisme del gradient per a funcions quadratiques   
%                                 requereix 'dds' previ amb "Q","c" i "eps"
% Inicialitzacions
kse=fopen('gradqm.res','w');
nit=300;                    % nombre maxim d'iteracions
n=length(c);
vaps=eig(Q);
a=vaps(1);                  % mes petit valor propi de Q
A=vaps(n);                  % mes gran valor propi de Q
fstc=((A-a)/(A+a))^2        % fita superior de la taxa de convergencia
fprintf(kse,'fstc=%g\n',fstc); 
xopt=-(Q\c);                % x* (punt solucio)
x0=zeros(n,1);
x=x0;
erran=.5*(x-xopt)'*Q*(x-xopt)
g=Q*x+c;
ngrq=g'*g;                  % norma del gradient al quadrat
d=-g;
%                            Iteracions ("nit" com a molt)
ite=0;
while ((ite<nit) & (ngrq>eps))
   ite=ite+1
   qd=Q*d;
   alf=ngrq/(d'*qd);
   x=x+alf*d;
   g=g+alf*qd;
   ngrq=g'*g
   d=-g;
%                             escriptura a arxiu
   fprintf(kse,' ite= %i ngrq= %g\n',ite,ngrq);
   fprintf(kse,' x:\n');
   for i=1:n
      fprintf(kse,' %g',x(i)); 
   end
   fprintf(kse,'\n');
end
fclose(kse);
