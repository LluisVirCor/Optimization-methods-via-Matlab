 % Dades problema estimacio estat (cas 45)
 % Numero de nusos
 nus = 6; 
 % Nus de referencia on s'ha pres la tensio
 nref = 6; 
 % Tensio mesurada al nus de referencia: e_nref=tref f_nref=0
 tref =   1.0005;
 % Dades sobre les linies: nlin, r, x, b 
 linia=[ 1  0.0011  0.0089  0.0244 ;        % linia  1
         2  0.0042  0.0206  0.0571 ;        % linia  2
         3  0.0019  0.0209  0.0519 ;        % linia  3
         4  0.0037  0.0183  0.0510 ;        % linia  4
         5  0.0021  0.0164  0.0427 ;        % linia  5
         6  0.0182  0.0885  0.2475 ;        % linia  6
         7  0.0118  0.0921  0.2330 ;        % linia  7
         8  0.0185  0.0890  0.0510 ];       % linia  8
 % Dades sobre les mesures de potencia: 
 % nus origen, nus desti, linia, pot. act., pot. react.
 mesura=[ 2  1  1 -0.2655  0.5855 ;
          2  3  6  0.1252 -0.0716 ;
          3  2  7 -0.1170 -0.1753 ;
          3  6  4 -0.6081 -0.6910 ;
          4  3  5  0.5294  0.4914 ;
          4  6  3 -0.0711 -0.2615 ;
          5  1  8 -0.2613  0.2375 ;
          5  6  2 -1.0581  0.4619 ;
          6  5  2  1.0627 -0.4964 ;
          6  3  4  0.6080  0.6536 ;
          6  4  3  0.0704  0.2068 ];
 % Punt inicial: x0=[1 1 1 1 1 1 0 0 0 0 0 0];
 x=[ones(nus,1);zeros(nus,1)];
