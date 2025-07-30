%
% Funcio que implementa l'exploracio lineal del l'algoritme A6.3.1mod del
% Dennis-Schnabel (comprova A.G.1 i A.G.2).
%   -> ENTRADA
%       xc - punt actual
%       fc - valor de la funcio objectiu
%       vfun - nom de la funcio objectiu
%       gc - gradient de la funcio objectiu
%       grad - nom de la funcio q calcula el gradient
%       p - direccio
%       Sx - Invers del valor t\ufffdpic de les variables
%       maxstep - longitud de passa maxima
%       steptol - tolerancia de la passa
%       beta - constant de la 2a cond. d'Armijo-Goldstein
%       varargin - parametres opcionals q utilitzen la fo i el gradient
%   <- SORTIDA
%       xplus - nou punt
%       fplus - nou valor de la funcio objectiu
%       gplus - valor actualitzat del gradient
%       retcode - (=0)-> punt acceptable (=1) -> punt no acceptable
%       maxtaken - (=1) -> passa maxima efectuada (si no =0)
%       l - longitud de passa optima
%       nf - num. avaluacions de la f.o.
%       ng - num. avaluacions del gradient
%

% Condicions d'Armijo-Goldstein:
% (AG1) f(x+lambda*d) <= f(x) + alfa*lambda*grad(x)*d     0 < alfa < 1
% (AG2) grad(x+lambda*d)*d >= beta*grad(x)*d           alfa < beta < 1

function [xplus,fplus,gplus,retcode,maxtaken,l,nf,ng]= ...
    lsearch(xc,fc,vfun,gc,grad,p,Sx,maxstep,steptol,beta,varargin)

nf = 0;
ng = 0;
maxtaken=0;
retcode=2;
alfa=10^(-4);
steplen=norm(diag(Sx)*p);
if (steplen > maxstep)
   p=p*(maxstep/steplen);
   steplen=maxstep;
end   
initslope=gc'*p;
%rellength=max(abs(p')/max(abs(xc'),(1/Sx')));
rellength = max(abs(p)./(max(xc,ones(length(Sx),1)./Sx)));
minl=steptol/rellength;
l=1;
cond1=(1==1);

while cond1
      xplus=xc+l*p;
      fplus=feval(vfun,xplus,varargin{:});
      nf = nf + 1;
      if (fplus <= fc+alfa*l*initslope)         % si AG1
         gplus=feval(grad,xplus,varargin{:});
         ng = ng + 1;
         newslope=gplus'*p;
         if (newslope < beta*initslope)         % no AG2
            if ((l == 1) & (steplen < maxstep))
               maxl=maxstep/steplen;
               fplus=feval(vfun,xplus,varargin{:});
               nf = nf + 1;
               cond2=(2==2);
               while cond2
                     lprev=l;
                     fplusprev=fplus;
                     l=min(2*l,maxl);
                     xplus=xc+l*p;
                     fplus=feval(vfun,xplus,varargin{:});
                     nf = nf + 1;
                     if (fplus <= fc+alfa*l*initslope)
                        gplus=feval(grad,xplus,varargin{:});
                        ng = ng + 1;
                        newslope=gplus'*p;
                     end    
                     if (fplus > fc +alfa*l*initslope) | (newslope >= beta*initslope) | (l >= maxl)
                        break
                     else
                        cond2=cond2;
                     end
               end
            end
            if (l < 1) | ((l > 1) & (fplus > fc+alfa*l*initslope))
               llo=min(l,lprev);
               ldiff=abs(lprev-l);
               if (l < lprev)
                  flo=fplus;
                  fhi=fplusprev;
               else
                  flo=fplusprev;
                  fhi=fplus;
               end   
               cond3=(3==3);
               while cond3 
                     lincr=((-1*newslope)*ldiff^2)/(2*(fhi-(flo+newslope*ldiff)));
                     if (lincr < 0.2*ldiff)
                        lincr=0.2*ldiff;
                     end
                     l=llo+lincr;
                     xplus=xc+l*p;
                     fplus=feval(vfun,xplus,varargin{:});
                     nf = nf + 1;
                     if (fplus > fc+alfa*l*initslope)
                        ldiff=lincr;
                        fhi=fplus;
                     else
                        gplus=feval(grad,xplus,varargin{:});
                        ng = ng + 1;
                        newslope=gplus'*p;
                        if (newslope < beta*initslope)
                           llo=l;
                           ldiff=ldiff-lincr;
                           flo=fplus;
                        end
                     end
                     if (newslope >= beta*initslope) | (ldiff < minl)
                        break
                     else
                        cond3=cond3;
                     end
               end
               if (newslope < beta*initslope)
                  fplus=flo;
                  xplus=xc+llo*p;
               end
            end   
         end      
         retcode=0;
         if (l*steplen > 0.99*maxstep)
            maxtaken=1;
         end
         return
      elseif (l < minl)
            retcode=1;
            gplus=feval(grad,xplus,varargin{:});
            ng = ng + 1;
            xplus=xc;
            return
      else
         if (l == 1) 
            ltemp=(-initslope)/(2*(fplus-fc-initslope));;
         else
            mat=[1/l^2 -1/lprev^2;-lprev/l^2 l/lprev^2];
            vec=[fplus-fc-l*initslope;fplusprev-fc-lprev*initslope];
            cons=1/(l-lprev);
            coefs=cons*mat*vec;
            a=coefs(1);
            b=coefs(2);
            disc=b^2-3*a*initslope;
            if (a == 0)
               ltemp=-initslope/(2*b);
            else
               ltemp=(-b+sqrt(disc))/(3*a);
            end
            if (ltemp > 0.5*l)
               ltemp=0.5*l;
            end
         end   
         lprev=l;
         fplusprev=fplus;
         if (ltemp <= 0.1*l)
            l=0.1*l;
         else
            l=ltemp;
         end
      end
      if (retcode < 2) 
         break
      else
         cond1=cond1;
      end
end
