[x1,y1]=textread('error_vs_t_btcs.txt','%f%f');
x11=x1./23;
y11=log(y1);
[x2,y2]=textread('ErrorVsTime_FTCS.txt','%f%f');
y22=log(y2);
[x3,y3]=textread('ErrorVsTime_ADI.txt','%f%f');
y33=log(y3);
plot(x11,y11,'r-',x2,y22,'b-',x3,y33,'g-')
xlabel('log(error)');
ylabel('log(error)');
xlabel('Time Iteration');
legend('BTCS','FTCS','ADI');
title('Convergence History for Stainless Steel')