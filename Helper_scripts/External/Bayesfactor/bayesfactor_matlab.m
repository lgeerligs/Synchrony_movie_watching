%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% MacLean et al. 2010 (correlation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 r = -0.3616346;
 n = 54;
 
 jzs_corbf(r,n) % BF10=3.85
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Kanai et al. in press (correlation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
r=0.4844199;
n=40;
 
jzs_corbf(r,n) % BF10=17.87
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Lleras et al. in press (partial correlation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
n=40;
p0=1;
p1=2;
r0=sqrt(0.6084);
r1=sqrt(0.6084408);
 
jzs_partcorbf(r0,r1,p0,p1,n) % BF10=0.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

