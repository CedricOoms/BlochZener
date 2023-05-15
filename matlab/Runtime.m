dt_vec=25:5:75;
runtime=[];
for dt=dt_vec
    runtime(end+1)=Qwalk_mod(dt);
end

figure();
plot(dt_vec,runtime);
xlabel("dt");
ylabel("runtime");
title("Runtime of the code as a function of the numerical time step size dt");