distance=sqrt((S(i).xd-(S(n+1).xd))^2 + (S(i).yd-(S(n+1).yd))^2);
dmax=sqrt((xm-(S(n+1).xd))^2 + (ym-(S(n+1).yd))^2);
cost = 1*(S(i).E/Eo)+ 1/7*distance/dmax;
if(temp_rand <=((p/(1-p*mod(r,round(1/p))))*cost))