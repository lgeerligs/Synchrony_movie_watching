function[bf10] = jzs_corbf(r0,r1,p0,p1,n)
bf10=quadgk(@(g)int(r1,p1,n,g),0,inf)/quadgk(@(g)int(r0,p0,n,g),0,inf);
