clear
clc
results=[];
k=5;
G=2000;
rep=3;
popN=40;
option2='soft';
%------------JD2001-------------
cd Liu_Tang_Yao
JD2001
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_2001_data.mat
instances=JD2001; 
N=size(instances,1); 
D=0.014*JDD2001;
vehicle_cap=capacity;
vehicle_num=K;
for I=1:3
    for n=1:rep
        ['JD_2001 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD2001_miu(I).data,JD2001_alpha(I).data,JD2001_beta(I).data,D);
    end
end
save MCEDA_DKE_JD2001_data.mat nd_pop nd_fitness
%------------JD2002-------------
cd Liu_Tang_Yao
JD2002
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_2002_data.mat  
instances=JD2002;
N=size(instances,1); 
D=0.014*JDD2002; 
vehicle_cap=capacity; 
vehicle_num=K; 
for I=1:3 
    for n=1:rep
        ['JD_2002 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD2002_miu(I).data,JD2002_alpha(I).data,JD2002_beta(I).data,D);
    end
end
save MCEDA_DKE_JD2002_data.mat nd_pop nd_fitness
%------------JD2003-------------
cd Liu_Tang_Yao
JD2003
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_2003_data.mat
instances=JD2003;
N=size(instances,1); 
D=0.014*JDD2003; 
vehicle_cap=capacity; 
vehicle_num=K; 
for I=1:3
    for n=1:rep
        ['JD_2003 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD2003_miu(I).data,JD2003_alpha(I).data,JD2003_beta(I).data,D);
    end
end
save MCEDA_DKE_JD2003_data.mat nd_pop nd_fitness
%------------JD2004-------------
cd Liu_Tang_Yao
JD2004
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load JD_2004_data.mat
instances=JD2004;
N=size(instances,1); 
D=0.014*JDD2004; 
vehicle_cap=capacity; 
vehicle_num=K;
for I=1:3
    for n=1:rep
        ['JD_2004 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,JD2004_miu(I).data,JD2004_alpha(I).data,JD2004_beta(I).data,D);
    end
end
save MCEDA_DKE_JD2004_data.mat nd_pop nd_fitness
