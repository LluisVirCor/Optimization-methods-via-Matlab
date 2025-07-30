%
% Calcul del gradient de la funcio d'estimacio d'estat 
% en el punt x
%
% Utilitza les funcions pot_act() i pot_rea()
%
function g=grad_estest(x,linia,mesura,nus,nref,tref)
[n,m]=size(mesura);

g=zeros(length(x),1);

for i=1:n
    ori=mesura(i,1);    % nus origen
    des=mesura(i,2);    % nus desti
    lin=mesura(i,3);    % linia
    
    % Recuperacio var: ek,fk i el,fl
    if (ori ~= nref)
        ek=x(ori);
        fk=x(nus+ori);
    else
        ek=tref;
        fk=0;
    end
    if (des ~= nref)
        el=x(des);
        fl=x(nus+des);
    else
        el=tref;
        fl=0;
    end
    
    % Calcul conductancia i susceptancia
    rkl=linia(lin,2);
    xkl=linia(lin,3);
    bkl=linia(lin,4);
    
    den = (rkl*rkl + xkl*xkl);
    ckl = rkl/den;
    skl = xkl/den;
    dkl = 2*skl - bkl;
    
    % Calcul pkl i qkl
    [pkl,aek,ael,afk,afl] = pot_act(ek,el,fk,fl,ckl,skl);
    [qkl,rek,rel,rfk,rfl] = pot_rea(ek,el,fk,fl,ckl,skl,dkl);
    
    % Valor del gradient
    if (ori ~= nref)
        g(ori) = g(ori) + ...
            2*((pkl-mesura(i,4))*aek + (qkl-mesura(i,5))*rek);
        g(nus+ori) = g(nus+ori) + ...
            2*((pkl-mesura(i,4))*afk + (qkl-mesura(i,5))*rfk);
    end
    if (des ~= nref)
        g(des) = g(des) + ...
            2*((pkl-mesura(i,4))*ael + (qkl-mesura(i,5))*rel);
        %CANVIS
        g(nus+des) = g(nus+des)+...
            2*((pkl-mesura(i,4))*ael+(qkl-mesura(i,5))*rel);
    end
end
