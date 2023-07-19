n_dt_vec=2:2:100;
runtime=[];
dt=[];
for n_dt=n_dt_vec
    [runtime_i,dt_i]=Qwalk_mod(n_dt);
    runtime(end+1)=runtime_i;
    dt(end+1)=dt_i;
end

figure();
plot(dt,runtime);
xlabel("dt");
ylabel("runtime");
title("Runtime of the code as a function of the numerical time step size dt");