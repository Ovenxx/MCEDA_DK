clear
clc
results=[];
k=5;
G=2000;
rep=3;
popN=40;
option2='soft';
%------------cdp101-------------
cd benchmark_Solomon
cdp101
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load cdp101_data.mat
instances=cdp101;
N=size(instances,1);
for I=1:5 
    vehicle_cap=capacity; 
    vehicle_num=K;
    for n=1:rep
        ['cdp101 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,cdp101_miu(I).data,cdp101_alpha(I).data,cdp101_beta(I).data);
    end
end
save MCEDA_DKE_cdp101_data.mat nd_pop nd_fitness
%------------cdp201-------------
cd benchmark_Solomon
cdp201
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load cdp201_data.mat
instances=cdp201;
N=size(instances,1); 
for I=1:5 
    vehicle_cap=capacity; 
    vehicle_num=K;
    for n=1:rep
        ['cdp201 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,cdp201_miu(I).data,cdp201_alpha(I).data,cdp201_beta(I).data);
    end
end
save MCEDA_DKE_cdp201_data.mat nd_pop nd_fitness 
%------------rdp101-------------
cd benchmark_Solomon
rdp101
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load rdp101_data.mat
instances=rdp101;
N=size(instances,1); 
for I=1:5 
    vehicle_cap=capacity; 
    vehicle_num=K;
    for n=1:rep
        ['rdp101 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,rdp101_miu(I).data,rdp101_alpha(I).data,rdp101_beta(I).data);
    end
end
save MCEDA_DKE_rdp101_data.mat nd_pop nd_fitness 
%------------rdp201-------------
cd benchmark_Solomon
rdp201
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load rdp201_data.mat
instances=rdp201;
N=size(instances,1); 
for I=1:5 
    vehicle_cap=capacity; 
    vehicle_num=K;
    for n=1:rep
        ['rdp201 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,rdp201_miu(I).data,rdp201_alpha(I).data,rdp201_beta(I).data);
    end
end
save MCEDA_DKE_rdp201_data.mat nd_pop nd_fitness 
%------------rcdp101-------------
cd benchmark_Solomon
rcdp101
current_dir=pwd;
parent_dir=fileparts(current_dir);
cd(parent_dir);
load rcdp101_data.mat
instances=rcdp101;
N=size(instances,1); 
for I=1:5 
    vehicle_cap=capacity; 
    vehicle_num=K;
    for n=1:rep
        ['rcdp101 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,rcdp101_miu(I).data,rcdp101_alpha(I).data,rcdp101_beta(I).data);
    end
end
save MCEDA_DKE_rcdp101_data.mat nd_pop nd_fitness 
%------------rcdp201-------------
cd benchmark_Solomon
rcdp201
current_dir=pwd;parent_dir=fileparts(current_dir);
cd(parent_dir);
load rcdp201_data.mat
instances=rcdp201;
N=size(instances,1); 
for I=1:5 
    vehicle_cap=capacity; 
    vehicle_num=K; 
    for n=1:rep
        ['rcdp201 ',num2str(I),num2str(n)]
        [nd_pop(I,n).data,nd_fitness(I,n).data]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,rcdp201_miu(I).data,rcdp201_alpha(I).data,rcdp201_beta(I).data);
    end
end
save MCEDA_DKE_rcdp201_data.mat nd_pop nd_fitness
