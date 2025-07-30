function f=estest(x,linia,mesura,nus,nref,tref)
% Calcula el valor de la funcio d'estimacio d'estat f donat un punt x
[n,m]=size(mesura);

f=0.0;
for i=1:n
    ori=mesura(i,1);    % nus origen
    des=mesura(i,2);    % nus desti
    lin=mesura(i,3);    % linia
    
    % Recuperacio ek,fk i el,fl del vector x()
    %CANVIS
  if (ori~=nref)
    	ek=x(ori); 
    	fk=x(ori+6);
	else
		ek=tref;
		fk=0;
	end
    
    if (des~=nref)
        el=x(des);
        fl=x(des+6);
    else
        el=tref;
        fl=0;
    end
    
    % Calcul conductancia i susceptancia
    rkl=linia(lin,2);
    xkl=linia(lin,3);
    bkl=linia(lin,4);

    den=(rkl*rkl + xkl*xkl);
    ckl=rkl/den;
    skl=xkl/den;
    
    % Calcul pkl i qkl
    %CANVIS
    aek= 2*ckl*ek - ckl*el - skl*fl;
    ael= -ckl*ek + skl*fk;
    afk= skl*el + 2*ckl*fk - ckl*fl; 
    afl= -skl*ek - ckl*fk;
    pkl=.5*(ek*aek+el*ael+fk*afk+fl*afl);
  
    dkl=2*skl-bkl;
    rek= dkl*ek - skl*el + ckl*fl;
    rel= -skl*ek -ckl*fk;
    rfk= -ckl*el + dkl*fk - skl*fl;
    rfl= ckl*ek - skl*fk;
    qkl= .5*(ek*rek + el*rel + fk*rfk + fl*rfl);
    
    % Valor de la funcio objectiu
    f= f +(pkl-mesura(i,4))^2 + (qkl-mesura(i,5))^2;
end
