%
% Minimitzacio d'una funcio qualsevol sense constriccions 
% Metode del GRADIENT CONJUGAT. 
% Utilitza la funcio lsearch.m d'exploracio lineal
%   -> ENTRADA
%       fun - Funcio a avaluar( f=fun(x,varargin))
%       grd - Gradient de fun ( g=grd(x,varargin))
%       x - Punt inicial
%       varargin - Arguments extres necessaris per avaluar la funcio
%   <- SORTIDA
%       x - Punt solucio
%       f - valor de la funcio objectiu
%       nite - num. iterarcions metode del gradient
%       nfun - num. crides a la funcio objectiu
%       ngrd - num. crides al calcul del gradient
%       xfob - valor de la funcio objectiu a cada iteracio
%       xnor - norma del gradient a cada iteracio
%
%   Crida:  [x,fval,nite,nfun,ngrd] = gradconj_fqsevol(@fun,@grd,x,varargin)
%           (per la funcio especifica d'estimacio d'estat)
%           [x,fval,nit,nfu,ngr] = ...
%              gradconj_fqsevol(@estest,@grad_estest,x0,linia,mesura,nus,nref,tref)
%
function [x,f,nite,nfun,ngrd,xfob,xnor] = ...
    gradconj_fqsevol(fun,grd,x,varargin)

x0=x(:);                    % Per treballar amb vectors columna
n=length(x0);               % Dimensio del problema

% Inicialitzacio parametres
mxite = 1000;               % Num. max. iteracions
mxfun = 1000;               % Num. max. avaluacio fun.obj.
mxgrd = 100;                % Num. max. avaluacio gradient

tol_opt = 1e-6;             % Tolerancia al punt optim

xfob = zeros(mxite,1);      % Valor de la fun. obj. a cada iteracio
xnor = zeros(mxite,1);      % Norma al quadrat del gradient a cada ite.

nite = 1;                   % Num. iteracions

uns=ones(n,1);              % Inversa del valor tipic de les variables
maxstep=(10^3)*max(norm(diag(uns)*x0),norm(diag(uns)));
maxstep=0.1;
maxstep
%%%%%%%%%%%%%%%%%%%%%   Algoritme   %%%%%%%%%%%%%%%%%%%%%
f = feval(fun,x,varargin{:});   % Avaluacio funcio objectiu
g = feval(grd,x,varargin{:});   % Avaluacio gradient
nfun = 1;                   % Num. avaluacions de la fun. obj.
ngrd = 1;                   % Num. avaluacions del gradient
norma = g'*g;
xfob(nite) = f;
xnor(nite) = norma;
disp(sprintf(' ite:    fun. obj.       norma   nf  ng      passa'));
disp(sprintf(' %3i: %12.5f   %9.3g  %3i %3i ',nite,f,norma,nfun,ngrd));
d = -g;
while ((nite<mxite) & (nfun<mxfun) & (ngrd<mxgrd) & (norma > tol_opt))
    nite = nite + 1;
    
    %CANVI
    normaux = g*g';

    % Passa: Exploracio lineal 
    [xplus,fplus,gplus,retcode,maxtaken,alf,nf,ng]= ...
        lsearch(x,f,fun,g,grd,d,uns,maxstep,eps^(2/3),0.9,varargin{:});
    nfun = nfun + nf;
    ngrd = ngrd + ng;
    % Criteri d'optimalitat
    den = norma;
    % Nou punt
    x(:) = xplus;
    % Nou gradient
    g = gplus;
    % Criteri d'optimalitat
    norma = g'*g;
    
    %CANVIS
    %Beta
    beta = g*g'/normaux
    
    % Direccio
    if rem(nite,n) == 0
       d=-g;
    else
       d = -g+beta*d
    end
    % Traça
    f=fplus;

    xfob(nite) = f;
    xnor(nite) = norma;

    disp(sprintf(' %3i: %12.5f   %9.3g  %3i %3i  %9.5f ',...
        nite,f,norma,nfun,ngrd,alf));
end

if (norma <= tol_opt)
    disp(sprintf('Solucio dins tolerancies'));
elseif (nite>=mxite)
    disp(sprintf('Num. max. iteracions exhaurit. Incrementa mxite.'));
elseif (nfun >= mxfun)
    disp(sprintf('Num. max. d''avaluacions de la funcio objectiu exhaurit.'));
    disp(sprintf('Incrementa mxfun.'));
elseif (ngrd >= mxgrd)
    disp(sprintf('Num. max. d''avaluacions del gradient exhaurit.'));
    disp(sprintf('Incrementa mxgrd.'));
end
