n_dt_vec=20:2:200;
runtime=[];
dt=[];
for n_dt=n_dt_vec
    runtime_dt=0;
    for i=1:5
        [runtime_i,dt_i]=Qwalk_mod(n_dt);
        runtime_dt=runtime_dt+runtime_i;
    end
    runtime(end+1)=runtime_dt/5;
    dt(end+1)=dt_i;
end

figure();
plot(dt,runtime);
xlabel("dt");
ylabel("runtime");
title("Runtime of the code as a function of the time step size dt");