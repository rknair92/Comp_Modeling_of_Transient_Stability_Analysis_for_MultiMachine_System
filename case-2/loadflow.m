function Pe=loadflow(t,tfi,tfc,Yred,Yredf,Yredpf,E,del);


if t < tfi
G=real(Yred);
B=imag(Yred);
elseif (t >=tfi) & (t < tfc) 
G=real(Yredf);
B=imag(Yredf);
else
G=real(Yredpf);
B=imag(Yredpf); 
end
 
n=length(G);
for p=1:n
    sum=0;
    for q=1:n
        if p~=q
        sum=sum+E(q)*(B(p,q)*sin(del(p)-del(q))+G(p,q)*cos(del(p)-del(q)));
        end
    end
    Pe(p)=sum*E(p)+E(p)*E(p)*G(p,p);
end

