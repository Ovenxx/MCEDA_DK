clear
clc
results=[];
k=5;
G=2000;
rep=50;
popN=40;
option2='soft';
%------------JD4001-------------
cd Liu_Tang_Yao
JD4001
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_4001_data.mat
instances=JD4001;
N=size(instances,1); 
D=0.014*JDD4001; 
vehicle_cap=capacity; 
vehicle_num=K; 
for I=1:3 
    for n=1:rep
        ['JD_4001 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD4001_miu(I).data,JD4001_alpha(I).data,JD4001_beta(I).data,D);
    end
end
save MCEDA_DKE_JD4001_data.mat nd_pop nd_fitness
%------------JD4002-------------
cd Liu_Tang_Yao
JD4002
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_4002_data.mat 
instances=JD4002;
N=size(instances,1); 
D=0.014*JDD4002; 
vehicle_cap=capacity; 
vehicle_num=K; 
for I=1:3
    for n=1:rep
        ['JD_4002 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD4002_miu(I).data,JD4002_alpha(I).data,JD4002_beta(I).data,D);
    end
end
save MCEDA_DKE_JD4002_data.mat nd_pop nd_fitness
%------------JD4003-------------
cd Liu_Tang_Yao
JD4003
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_4003_data.mat
instances=JD4003;
N=size(instances,1);
D=0.014*JDD4003;
vehicle_cap=capacity; 
vehicle_num=K;
for I=1:3
    for n=1:rep
        ['JD_4003 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD4003_miu(I).data,JD4003_alpha(I).data,JD4003_beta(I).data,D);
    end
end
save MCEDA_DKE_JD4003_data.mat nd_pop nd_fitness
%------------JD4004-------------
cd Liu_Tang_Yao
JD4004
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_4004_data.mat
instances=JD4004;
N=size(instances,1);
D=0.014*JDD4004;
vehicle_cap=capacity; 
vehicle_num=K; 
for I=1:3 
    for n=1:rep
        ['JD_4004 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD4004_miu(I).data,JD4004_alpha(I).data,JD4004_beta(I).data,D);
    end
end
save MCEDA_DKE_JD4004_data.mat nd_pop nd_fitness
