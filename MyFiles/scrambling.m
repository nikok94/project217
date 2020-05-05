function scramb_seq = scrambling(c_init,E)
scramb_seq=zeros(1,E);
c_init=uint32(c_init);
Nc=1600;
x1=zeros(1,Nc+E);
x1(1)=1;
x2=zeros(1,Nc+E);
for i=0:30
    x2(i+1)=bitand(c_init,uint32(1));
    c_init=bitshift(c_init,-1);
end
for n=0:Nc+E-31-1
    x1(n+31+1)=mod(x1(n+3+1)+x1(n+1),2);
    x2(n+31+1)=mod(x2(n+3+1)+x2(n+2+1)+x2(n+1+1)+x2(n+1),2);
end

for n=0:E-1
    scramb_seq(n+1)=mod(x1(n+Nc+1)+x2(n+Nc+1),2);
end