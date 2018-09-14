function[bf10] = jzs_corbf(r,n)
p=1;
bf10=sqrt((n/2))/gamma(1/2)*quadgk(@(g)int(r,p,n,g),0,inf);
