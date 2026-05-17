%%
clear all
clc
instance={'cdp101','cdp201','rdp101','rdp201','rcdp101','rcdp201'};
for x=1:6
    load(string(strcat('MCEDA_DK_data\MCEDA_DK_',instance(x),'_data.mat')));
    % load(string(strcat('MCEDA_DKE_',instance(x),'_data.mat')));
    MCEDA_DKE_fitness=nd_fitness;
    load(string(strcat('MD_MOLS_data\MD_MOLS_',instance(x),'_data.mat')));
    MD_DKE_MOLS=nd_fitness;
    load(string(strcat('MD_MOMA_data\MD_MOMA_',instance(x),'_data.mat')));
    MD_DKE_MOMA=nd_fitness;
    for I=1:5
        s=[];
        N1 = length(MCEDA_DKE_fitness(I,:));
        for n=1:N1
            MCEDA_DKE_data(n).PF=MCEDA_DKE_fitness(I,n).data(2000).PF;
            s=[s;MCEDA_DKE_data(n).PF];
        end
        N2 = length(MD_DKE_MOLS(I,:));
        for n=1:N2
            MD_MOLS_data(n).PF=MD_DKE_MOLS(I,n).data(2000).PF;
            s=[s;MD_MOLS_data(n).PF];
        end
        N3 = length(MD_DKE_MOMA(I,:));
        for n=1:N3
            MD_MOMA_data(n).PF=MD_DKE_MOMA(I,n).data(2000).PF;
            s=[s;MD_MOMA_data(n).PF];
        end
        r=[1.5*max(s(:,1)),1.5*max(s(:,2))];
        for n=1:N1
            hv1(I,n)=hypervolume(MCEDA_DKE_data(n).PF,r,10^5);
        end
        for n=1:N2
            hv2(I,n)=hypervolume(MD_MOLS_data(n).PF,r,10^5);
        end
        for n=1:N3
            hv3(I,n)=hypervolume(MD_MOMA_data(n).PF,r,10^5);
        end
        HV1(x,I)=mean(hv1(I,:));
        HV2(x,I)=mean(hv2(I,:));
        HV3(x,I)=mean(hv3(I,:));
    end
end
HV1_column=HV1';
HV1_column=HV1_column(:);
HV2_column=HV2';
HV2_column=HV2_column(:);
HV3_column=HV3';
HV3_column=HV3_column(:);
rHV=[HV1_column,HV2_column,HV3_column];