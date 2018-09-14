        function[t] = int(r,p,n,g)
    t=exp( ((n-p-1)./2).*log(1+g) + (-(n-1)./2).*log(1+(1-r.^2).*g) +    (-3./2)*log(g) -n./(2.*g));