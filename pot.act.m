pot.act.m


function[pact,aek,ael,afk,afl]=...
pot_act(etk,etl,ftk,ftl,cil,sil)

M=[2+cil -cil   0    -sil;
   -cil   0    sil    0;
    0     sil  2*cil -cil;
   -sil   0    -cil   0];

v=[etk;etl;ftk;ftl];
aek=M(1,:)*v;
ael=M(2,:)*v;
afk=M(3,:)*v;
afl=M(4,:)*v;
pact=0.5*v'*M*v;
