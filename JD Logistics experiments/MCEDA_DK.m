% Multi-Consensus Based Estimation of Distribution Algorithm with Domain Knowledge for Multi-Objective Green Vehicle Routing
function [nd_pop,nd_fitness,time]=MCEDA_DK(instances,vehicle_cap,vehicle_num,N,k,G,popN,Miu,ALpha,BEta,D)
tic
lambda=0.90;
customer_dem=instances(:,2:3); 
customer_dem(:,1)=-customer_dem(:,1); 
customer_time=instances(:,4:6); 
scale_m=ceil(sqrt(N-1)); 
scale_delta=ceil(sqrt(N-1)/3); 
neisize=scale_m-scale_delta:scale_m+scale_delta;
% Initialization Stage
SC=proba_transform(D,N,lambda);
TC=time_compatibility(customer_time,D,N); 
H=SC.*TC; 
flag=0;
[Population,Fitness]=construct_population(customer_time,D,customer_dem,vehicle_cap,vehicle_num,H,N,k*popN,Miu,ALpha,BEta,flag); 
[CM1,CM2]=MultiConsensus_Computation(Population,Fitness,k*popN,N);
[Population,Fitness]=survival_of_the_fittest(Population,Fitness,popN);
nd_fitness=struct('PF', cell(1, G));
g_counter=0;
% Evolution Stage
for g=1:G

    [new_Population,new_Fitness]=MCE(customer_time,D,customer_dem,vehicle_cap,vehicle_num,H,CM1,CM2,N,popN,Miu,ALpha,BEta,flag);

    [new_Population_PC_LNS,new_Fitness_PC_LNS]=PC_LNS(new_Population,new_Fitness,customer_time,D,N,H,CM1,CM2,neisize,customer_dem,vehicle_cap,vehicle_num,popN,Miu,ALpha,BEta);
    [Population_PC_LNS,Fitness_PC_LNS]=PC_LNS(Population,Fitness,customer_time,D,N,H,CM1,CM2,neisize,customer_dem,vehicle_cap,vehicle_num,popN,Miu,ALpha,BEta);
    
    Population=[new_Population;Population;new_Population_PC_LNS;Population_PC_LNS];
    Fitness=[new_Fitness;Fitness;new_Fitness_PC_LNS;Fitness_PC_LNS];
    
    [Population,Fitness]=survival_of_the_fittest(Population,Fitness,popN);
    
    [Population_TFR,Fitness_TFR]=TFR(Population,Fitness,customer_dem,customer_time,vehicle_cap,D,popN,N,vehicle_num,Miu,ALpha,BEta);

    Population=[Population_TFR;Population];
    Fitness=[Fitness_TFR;Fitness];
    
    [Population,Fitness]=survival_of_the_fittest(Population,Fitness,popN);
    
    [CM1,CM2]=MultiConsensus_Computation(Population,Fitness,popN,N);
    [nd_pop,nd_fitness(g).PF]=nondominated_output(Population,Fitness);
    
    if min(nd_fitness(g).PF(:,1))>10^8 
        g_counter=g_counter+1;
    elseif flag==1
        flag=0;
    end 
    if g_counter>100
        flag=1;
        g_counter=0;
    end

end
time=toc;
end

% Local Functions
function [nd_pop,nd_fitness]=nondominated_output(Population,Fitness)
[sorted_rank,sorted_pop,sorted_fitness]=non_dominated_sorting(Population,Fitness);
for i=1
    POP=sorted_pop(sorted_rank==i,:);
    FIT=sorted_fitness(sorted_rank==i,:);
    crowding=crowding_distance(FIT);
    [~,sorted_indices]=sort(crowding,'descend');
    nd_pop=POP(sorted_indices,:);
    nd_fitness=FIT(sorted_indices,:);
end
end

function y=Nonlinear_modification(x,gre)
y=(1-x.^(1/gre)).^gre;
end

function H=proba_transform(D,N,lambda)
epislon=10^(-10);
for n1=1:N
    if n1==1
        Dslace=D(n1,n1+1:end);
    elseif n1==N
        Dslace=D(n1,1:n1-1);
    else
        Dslace=[D(n1,1:n1-1),D(n1,n1+1:end)];
    end
    [~,R1]=sort(Dslace);
    R=[];
    for nnn=1:N-1
            R(nnn)=find(R1==nnn);
    end
    if n1==1
        R=[0,R];
    elseif n1==N
        R=[R,0];
    else
        R=[R(1:n1-1),0,R(n1:end)];
    end
    maxDslace=max(Dslace);
    for n2=1:N
        if n2~=n1
            H(n1,n2)=lambda^(R(n2)-1)*(maxDslace-D(n1,n2)+epislon);
        else
            H(n1,n2)=0;
        end
    end
    H(n1,:)=H(n1,:)/sum(H(n1,:)); 
end
end

function TC=time_compatibility(customer_time,D,N)
Vmax=1000; 
Vmin=200; 
delta1=10^(-1); 
delta2=0;
epsilon=10^(-20);
ear=customer_time(:,1);
lat=customer_time(:,2);
ser=customer_time(:,3);
TC=zeros(N,N);
for i=1:N
    if i==1 
        for j=i+1:N
            TC(i,j)=(max(min(D(i,j)/Vmin,lat(j)),ear(j))-ear(j))/(lat(j)-ear(j))-delta1*(ear(j)-D(i,j)/Vmin);
        end
        TC(1,2:end)=(TC(1,2:end)-min(TC(1,2:end))+epsilon)/(max(TC(1,2:end))-min(TC(1,2:end))+epsilon);
    else 
        for j=1:N
            if j~=i
                if ear(i)+ser(i)+D(i,j)/Vmax>lat(j)
                    TC(i,j)=delta2;
                elseif lat(i)+ser(i)+D(i,j)/Vmin>=ear(j)
                    TC(i,j)=1;
                else
                    TC(i,j)=1/(1+ear(j)-(lat(i)+ser(i)+D(i,j)/Vmin));
                end
            end
        end
    end
    TC(i,:)=TC(i,:)/sum(TC(i,:));
end
end

function i=RouletteWheelSelection(P)
P=P/sum(P);
C=cumsum(P);
i=sum(rand>C)+1;
end

function [Population,Fitness]=MCE(customer_time,D,customer_dem,vehicle_cap,vehicle_num,H,CM1,CM2,N,popN,Miu,ALpha,BEta,flag)
w=linspace(0,1,popN);
Population=zeros(popN,N+vehicle_num);
Fitness=zeros(popN,2);
for n=1:popN
    P=H.*(w(n)*CM1+(1-w(n))*CM2);
    [Population(n,:),Fitness(n,:)]=construct_population(customer_time,D,customer_dem,vehicle_cap,vehicle_num,P,N,1,Miu,ALpha,BEta,flag);
end
end

function [Population,Fitness]=construct_population(customer_time,D,customer_dem,vehicle_cap,vehicle_num,P,N,popN,Miu,ALpha,BEta,flag)
ear=customer_time(:,1);
lat=customer_time(:,2);
ser=customer_time(:,3);
Vmax=1200;
Tmax=max(lat);
dem_delivery=-customer_dem(:,1);
dem_pickup=customer_dem(:,2);
P(:,1)=0;
P_copy=P;
Population=zeros(popN,N+vehicle_num);
Fitness=zeros(popN,2);
for n=1:popN
    P=P_copy;
    i=1; 
    S=zeros(1,N+vehicle_num); 
    S(i)=1; 
    used_vehicle_num=1;
    arrt=0;
    Nd=0;
    Np=0;
    nd=0;
    np=0;
    present=false(1,N);
    present(S(S~=0))=true;
    set_remain=2:N;
    while ~all(present) 
        if flag==1 
            f=ALpha(S(i),set_remain).*exp(-(arrt+ser(S(i))-Miu(S(i),set_remain)).^2./(Tmax*BEta(S(i),set_remain)).^2);
            V=(1-f)*Vmax;
            judge=arrt+ser(S(i))+D(S(i),set_remain)./V<=lat(set_remain)';
            P(S(i),set_remain(judge==0))=0;
        end 
        j=RouletteWheelSelection(P(S(i),:));
        if  j==1
            i=i+1;
            S(i)=1;
            used_vehicle_num=used_vehicle_num+1;
            arrt=0;
            Nd=0;
            Np=0;
            nd=0;
            np=0;
            if used_vehicle_num>vehicle_num
                P=P_copy;
                i=1;
                S=zeros(1,N+vehicle_num);
                S(i)=1;
                used_vehicle_num=1;
                set_remain=2:N;
            end
            present=false(1,N);
            present(S(S~=0))=true;
        else
            feas1=(sum(dem_delivery(j)<=vehicle_cap-Nd+nd-np)==length(nd))&&(dem_pickup(j)<=vehicle_cap-Np);
            if feas1 
                i=i+1;
                S(i)=j;
                P(:,j)=0;
                Nd=Nd+dem_delivery(j);
                Np=Np+dem_pickup(j);
                nd=[nd,Nd];
                np=[np,Np];
                set_remain(set_remain==j)=[];
                if flag==1
                    ff=ALpha(S(i-1),j)*exp(-(arrt+ser(S(i-1))-Miu(S(i-1),j))^2/(Tmax*BEta(S(i-1),j))^2);
                    V=(1-ff)*Vmax;
                    arrt=max(ear(j),arrt+ser(S(i-1))+D(S(i-1),j)/V);
                end
            else 
                i=i+1;
                S(i)=1;
                used_vehicle_num=used_vehicle_num+1; 
                arrt=0;
                Nd=0;
                Np=0;
                nd=0;
                np=0;
                if used_vehicle_num>vehicle_num
                    P=P_copy;
                    i=1; 
                    S=zeros(1,N+vehicle_num); 
                    S(i)=1; 
                    used_vehicle_num=1;
                    set_remain=2:N;
                end
            end
            present=false(1,N);
            present(S(S~=0))=true;
        end
    end
    S(i+1)=1;
    Population(n,:)=S;
    Fitness(n,:)=Fitness_evaluation(Population(n,:),1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
end
end

function [Fitness,Penalty_ltw,Penalty_rtw]=Fitness_evaluation(Population,popN,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta)
lambda=10^8;
vdc=2000;
tdc1=10;
tdc2=10;
ear=customer_time(:,1);
lat=customer_time(:,2);
ser=customer_time(:,3);
Tmax=max(lat);
Penalty_ltw=zeros(1,popN);
Penalty_rtw=zeros(1,popN);
Vmax=1200;
Wv=3850;
lam=1/(44*737); k=0.2; Ne=38.33; Ed=4.70; gamma=1/(1000*0.4*0.9); bet=0.5*0.7*3.912*1.2041;
alp=9.81*(sin(0)+0.01*cos(0));
fi_e=lam*k*Ne*Ed;
fi_s=lam*gamma*bet;
fi_w=lam*gamma*alp;
rou=2.68; 
Fitness=zeros(popN,2);
for n1=1:popN
    s=Population(n1,:);
    s(s==0)=[];
    depot_index=find(s==1);
    ns=length(s);
    DIS=0; 
    CBE=0; 
    arrt=0;
    depot_i=0;
    for n2=1:ns-1
        ve=1; 
        if s(n2)==1
            Total_Duration_single(ve)=0;
            depot_i=depot_i+1;
            F_Loads=sum(dem_delivery(s(depot_index(depot_i):depot_index(depot_i+1))));
        else
            F_Loads=F_Loads-dem_delivery(s(n2))+dem_pickup(s(n2));
        end
        f=ALpha(s(n2),s(n2+1))*exp(-(arrt+ser(s(n2))-Miu(s(n2),s(n2+1)))^2/(Tmax*BEta(s(n2),s(n2+1)))^2); 
        V=(1-f)*Vmax; 
        DIS=DIS+D(s(n2),s(n2+1));
        tij=D(s(n2),s(n2+1))/V; 
        Total_Duration_single(ve)=Total_Duration_single(ve)+tij;
        if s(n2+1)==1
            arrt=0;
        else
            arrt=arrt+ser(s(n2))+tij;
        end
        Penalty_ltw(n1)=Penalty_ltw(n1)+max(ear(s(n2+1))-arrt,0); 
        arrt=max(ear(s(n2+1)),arrt); 
        Penalty_rtw(n1)=Penalty_rtw(n1)+max(arrt-lat(s(n2+1)),0); 
        arrt=min(arrt,lat(s(n2+1)));
        V_std=V/60;
        e_012=fi_e/V_std+fi_s*V_std^2+fi_w*(Wv+F_Loads*10e3); 
        CBE=CBE+e_012*D(s(n2),s(n2+1))*rou;
    end
    Fitness(n1,1)=DIS+vdc*(sum(s==1)-1)+tdc1*(sum(Total_Duration_single)+tdc2*max(Total_Duration_single));
    Fitness(n1,2)=CBE;
    Fitness(n1,:)=Fitness(n1,:)+lambda*Penalty_rtw(n1);
end
end

function dom=check_domination(a,b)
if all(a <= b) && any(a < b)
    dom = 1;
elseif all(b <= a) && any(b < a)
    dom = -1;
else
    dom = 0;
end
end

function [sorted_pop,sorted_fitness]=population_update(Population,Fitness)
[sorted_rank,sorted_pop,sorted_fitness]=non_dominated_sorting(Population,Fitness);
front_num=unique(sorted_rank);
for i=front_num'
    POP=[];
    FIT=[];
    POP=sorted_pop(sorted_rank==i,:);
    FIT=sorted_fitness(sorted_rank==i,:);
    crowding=crowding_distance(FIT);
    [~,sorted_indices]=sort(crowding,'descend');
    sorted_pop(sorted_rank==i,:)=POP(sorted_indices,:);
    sorted_fitness(sorted_rank==i,:)=FIT(sorted_indices,:);
end
end

function [sorted_rank,sorted_pop,sorted_fitness]=non_dominated_sorting(pop,fitness)
[M, ~] = size(pop);
n_dominated = zeros(M, 1);
dominating = cell(M, 1);
for i = 1:M-1
    for j = i+1:M
        dom=check_domination(fitness(i, :), fitness(j, :));
        if dom==1
            dominating{i} = [dominating{i}, j];
            n_dominated(j) = n_dominated(j) + 1;
        elseif dom==-1
            dominating{j} = [dominating{j}, i];
            n_dominated(i) = n_dominated(i) + 1;
        end
    end
end
rank = zeros(M, 1);
front = find(n_dominated == 0);
current_rank = 1;
while ~isempty(front)
    rank(front) = current_rank;
    next_front = [];
    for i = 1:length(front)
        idx = front(i);
        for j = 1:length(dominating{idx})
            k = dominating{idx}(j);
            n_dominated(k) = n_dominated(k) - 1;
            if n_dominated(k) == 0
                next_front = [next_front, k];
            end
        end
    end
    current_rank = current_rank + 1;
    front = next_front;
end
[sorted_rank, sorted_indices] = sort(rank);
sorted_pop = pop(sorted_indices, :);
sorted_fitness = fitness(sorted_indices, :);
end

function crowding=crowding_distance(fitness_front)
[N, M] = size(fitness_front);
crowding = zeros(N, 1);
for m = 1:M
    f = fitness_front(:, m);
    [sorted_f, sorted_indices] = sort(f);
    f_min = sorted_f(1);
    f_max = sorted_f(end);
    if f_max == f_min
        continue;
    end
    delta = zeros(N, 1);
    delta(sorted_indices(1)) = Inf;
    delta(sorted_indices(end)) = Inf;
    for i = 2:length(sorted_indices)-1
        original_idx = sorted_indices(i);
        delta(original_idx) = (sorted_f(i+1) - sorted_f(i-1)) / (f_max - f_min);
    end
    crowding = crowding + delta;
end
end

function [CM1,CM2]=MultiConsensus_Computation(Population,Fitness,popN,N)
epislon=10^(-20);
alpha=0.1;
CM1=zeros(N,N);
CM2=zeros(N,N);
fitness=Fitness(:,1);
fitness=(max(fitness)-fitness+epislon)/(max(fitness)-min(fitness)+epislon);
fitness=exp(fitness);
f1=fitness/sum(fitness);
fitness=Fitness(:,2);
fitness=(max(fitness)-fitness+epislon)/(max(fitness)-min(fitness)+epislon);
fitness=exp(fitness);
f2=fitness/sum(fitness);
for n=1:popN
    A=zeros(N,N);
    s=Population(n,:);
    s(s==0)=[];
    solutionL=length(s);
    for n1=1:solutionL-1
        A(s(n1),s(n1+1))=1;
    end
    CM1=CM1+f1(n)*A;
    CM2=CM2+f2(n)*A;
end
theta=zeros(N,N);
theta(CM1==0)=1;
theta=theta-eye(N);
CMc=CM1;
CMc(CMc==0)=inf;
Mmin=mean(min(CMc));
CM1=CM1+alpha*Mmin*theta;
for n=1:N
    CM1(n,:)=CM1(n,:)/sum(CM1(n,:));
end
theta=zeros(N,N);
theta(CM2==0)=1;
theta=theta-eye(N);
CMc=CM2;
CMc(CMc==0)=inf;
Mmin=mean(min(CMc));
CM2=CM2+alpha*Mmin*theta;
for n=1:N
    CM2(n,:)=CM2(n,:)/sum(CM2(n,:));
end
end

function [SLNS_Population,SLNS_Fitness,R]=SLNS(Population,Fitness,customer_time,D,popN,N,H,CM,neisize,customer_dem,vehicle_cap,vehicle_num,Miu,ALpha,BEta)
SLNS_Population=[];
SLNS_Fitness=[];
dem_delivery=-customer_dem(:,1); 
dem_pickup=customer_dem(:,2);
maxTry=5;
P=H.*CM; 
for n=1:popN
    individual=Population(n,:);
    individual(individual==0)=[];
    neinum=neisize(RouletteWheelSelection(1./neisize));
    neighbors=randperm(N-1,neinum)+1;
    destory_individual=individual;
    for n_nei=1:neinum
        destory_individual(destory_individual==neighbors(n_nei))=[];
    end
    for a=1:length(destory_individual)-1
        if destory_individual(a)==1&&destory_individual(a+1)==1
            destory_individual(a)=0;
        end
    end
    destory_individual(destory_individual==0)=[];
    for nnn=1:maxTry
        copy_neighbors=neighbors(randperm(length(neighbors)));
        copy_destory_individual=destory_individual;
        while ~isempty(copy_neighbors)
            [copy_destory_individual,copy_neighbors,F]=Repair(copy_destory_individual,copy_neighbors,P,customer_dem,vehicle_cap,2);
            if F==1
                break;
            end
        end
        if F==0
            break;
        end
    end
    if F==0
        R=0;
        repair_individual=copy_destory_individual;
        repair_fitness=Fitness_evaluation(repair_individual,1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
        sss=[repair_individual,zeros(1,N+vehicle_num-length(repair_individual))];
        SLNS_Population=[SLNS_Population;sss];
        SLNS_Fitness=[SLNS_Fitness;repair_fitness];
    else
        R=1;
        SLNS_Population=[SLNS_Population;Population(n,:)];
        SLNS_Fitness=[SLNS_Fitness;Fitness(n,:)];
    end
end
end

function [SLNS_Population,SLNS_Fitness,Population,Fitness,F,SUCC]=SLNS_route_descent(Population,Fitness,customer_time,D,popN,N,H,CM,customer_dem,vehicle_cap,vehicle_num,Miu,ALpha,BEta)
SLNS_Population=[];
SLNS_Fitness=[];
dem_delivery=-customer_dem(:,1);
dem_pickup=customer_dem(:,2);
maxTry=5;
P=H.*CM;
SUCC=0;
for n=1:popN
    individual=Population(n,:);
    individual(individual==0)=[];
    depot_index=find(individual==1);
    NV=length(depot_index)-1;
    tour_cap=zeros(1,NV);
    for ii=1:NV
        tour_cap(ii)=depot_index(ii+1)-depot_index(ii)-1;
    end
    [~,neighbors_depot]=min(tour_cap);
    ss=depot_index(neighbors_depot)+1;
    ee=depot_index(neighbors_depot+1)-1;
    neighbors=individual(ss:ee);
    destory_individual=individual;
    destory_individual(ss:ee+1)=[]; 
    for nnn=1:maxTry
        copy_neighbors=neighbors(randperm(length(neighbors)));
        copy_destory_individual=destory_individual;
        while ~isempty(copy_neighbors)
            [copy_destory_individual,copy_neighbors,F]=Repair(copy_destory_individual,copy_neighbors,P,customer_dem,vehicle_cap,2);
            if F==1
                break;
            end
        end
        if F==0
            break;
        end
    end
    if F==0
        repair_individual=copy_destory_individual;
        [repair_fitness,~,new_Penalty_rtw]=Fitness_evaluation(repair_individual,1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
        sss=[repair_individual,zeros(1,N+vehicle_num-length(repair_individual))];
        SLNS_Population=[SLNS_Population;sss];
        SLNS_Fitness=[SLNS_Fitness;repair_fitness];
        if new_Penalty_rtw==0
            SUCC=1;
            Population(n,:)=[repair_individual,zeros(1,N+vehicle_num-length(repair_individual))];
            Fitness(n,:)=repair_fitness;
        end
    else
        SLNS_Population=[SLNS_Population;Population(n,:)];
        SLNS_Fitness=[SLNS_Fitness;Fitness(n,:)];
    end
end
end

function [Population,Fitness]=PC_LNS(Population,Fitness,customer_time,D,N,H,CM1,CM2,neisize,customer_dem,vehicle_cap,vehicle_num,popN,Miu,ALpha,BEta)
eArchive_P=[];
eArchive_F=[];
dem_delivery=-customer_dem(:,1);
dem_pickup=customer_dem(:,2);
REP_index=[];
for n=1:popN
    r=rand;
    CM=r*CM1+(1-r)*CM2;
    [~,~,penalty_rtw_individual]=Fitness_evaluation(Population(n,:),1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta); 
    if penalty_rtw_individual>0||sum(Population(n,:)==1)==2 
        [Population(n,:),Fitness(n,:),R]=SLNS(Population(n,:),Fitness(n,:),customer_time,D,1,N,H,CM,neisize,customer_dem,vehicle_cap,vehicle_num,Miu,ALpha,BEta);
        if R==1 
            REP_index=[REP_index,n];
        end
    else 
        [Population_RD,Fitness_RD,Population(n,:),Fitness(n,:),F,SUCC]=SLNS_route_descent(Population(n,:),Fitness(n,:),customer_time,D,1,N,H,CM,customer_dem,vehicle_cap,vehicle_num,Miu,ALpha,BEta);
        if F==0&&SUCC==1 
            
        elseif F==1 
            [Population(n,:),Fitness(n,:),R]=SLNS(Population(n,:),Fitness(n,:),customer_time,D,1,N,H,CM,neisize,customer_dem,vehicle_cap,vehicle_num,Miu,ALpha,BEta);
            if R==1 
                REP_index=[REP_index,n];
            end
        elseif F==0&&SUCC==0 
            [Population(n,:),Fitness(n,:)]=SLNS(Population_RD,Fitness_RD,customer_time,D,1,N,H,CM,neisize,customer_dem,vehicle_cap,vehicle_num,Miu,ALpha,BEta);
            dom=check_domination(Fitness(n,:),Fitness_RD);
            if  dom==-1 
                Population(n,:)=Population_RD;
                Fitness(n,:)=Fitness_RD;
            elseif dom==0 
                eArchive_P=[eArchive_P;Population_RD];
                eArchive_F=[eArchive_F;Fitness_RD];
            end
        end
    end
end
Population(REP_index,:)=[];
Fitness(REP_index,:)=[];
Population=[Population;eArchive_P];
Fitness=[Fitness;eArchive_F];
end

function [destory_individual,neighbors,F]=Repair(destory_individual,neighbors,L,customer_dem,vehicle_cap,gre)
epislon=10^(-20);
dem_delivery=customer_dem(:,1);
dem_pickup=customer_dem(:,2);
F=0;
midnei=neighbors(1); 
neighbors(1)=[]; 
len=length(destory_individual);
p=zeros(1,len-1); 
for i=1:len-1
    p(i)=L(destory_individual(i),midnei)*L(midnei,destory_individual(i+1)); 
end
p=(1-((max(p)-p+epislon)./(max(p)-min(p)+epislon)).^(1/gre)).^gre;
T=0;
while (T==0)&&(sum(p)~=0) 
    insert_position=RouletteWheelSelection(p); 
    p(insert_position)=0;
    middi=destory_individual(1:insert_position);
    vnumber=sum(middi(middi==1)); 
    individual_copy=[destory_individual(1:insert_position),midnei,destory_individual(insert_position+1:end)];
    indexx=find(individual_copy==1);
    s=individual_copy(indexx(vnumber):indexx(vnumber+1)); 
    [flag,rofinf]=Feasibility_check(s,dem_delivery,dem_pickup,vehicle_cap); 
    if flag==1
        T=1; 
    else
        if rofinf==1 
            p(indexx(vnumber):indexx(vnumber+1)-1)=0;
        end
    end
end
if T==0
    F=1;
else
    destory_individual=[destory_individual(1:insert_position),midnei,destory_individual(insert_position+1:end)];
end
end

function [Population,Fitness]=TFR(Population,Fitness,customer_dem,customer_time,vehicle_cap,D,popN,N,vehicle_num,Miu,ALpha,BEta)
dem_delivery=customer_dem(:,1); 
ne_dem_delivery=-dem_delivery; 
dem_pickup=customer_dem(:,2);
max_fragment_length=3;
remove_n=[];
gre=2;
for n=1:popN
    ori_p=Population(n,:);
    if sum(Population(n,:)==1)>=3 
        [Population(n,:),Fitness(n,:)]=inter_TFR_exchange_SWAPstar(Population(n,:),Fitness(n,:),ne_dem_delivery,dem_pickup,customer_time,vehicle_cap,D,N,vehicle_num,gre,Miu,ALpha,BEta,max_fragment_length);
        [Population(n,:),Fitness(n,:)]=inter_TFR(Population(n,:),Fitness(n,:),ne_dem_delivery,dem_pickup,customer_time,vehicle_cap,D,N,vehicle_num,gre,Miu,ALpha,BEta,max_fragment_length);
    end
    [Population(n,:),Fitness(n,:)]=intra_TFR_exchange(Population(n,:),Fitness(n,:),dem_delivery,dem_pickup,customer_time,vehicle_cap,D,gre,Miu,ALpha,BEta,max_fragment_length);
    [Population(n,:),Fitness(n,:)]=intra_TFR(Population(n,:),Fitness(n,:),ne_dem_delivery,dem_pickup,customer_time,vehicle_cap,D,gre,Miu,ALpha,BEta,max_fragment_length);
    if sum(ori_p~=Population(n,:))==0
        remove_n=[remove_n,n];
    end
end
Population(remove_n,:)=[]; 
Fitness(remove_n,:)=[];
end

function [individual,fitness]=inter_TFR(individual,fitness,dem_delivery,dem_pickup,customer_time,vehicle_cap,D,N,vehicle_num,gre,Miu,ALpha,BEta,max_fragment_length)
s=individual;
index_of_depots=find(s==1);
number_of_depots=length(index_of_depots)-1;
tour_index1=randi(number_of_depots);
index_of_depots_t1=index_of_depots(tour_index1);
tour=s(index_of_depots_t1:index_of_depots(tour_index1+1));
number_of_customers=length(tour)-2;
max_fragment_length=min(number_of_customers,max_fragment_length);
fragment_length=randi(max_fragment_length); 
start_p=randi(number_of_customers-fragment_length+1)+1;
end_p=start_p+fragment_length-1;
fragment1=tour(start_p:end_p); 
xxx=tour;
xxx(start_p:end_p)=[];
cluster_d=zeros(1,number_of_depots);
for n=1:number_of_depots
    if n~=tour_index1
        midd=s(index_of_depots(n)+1:index_of_depots(n+1)-1);
        dd=0;
        for ii=1:number_of_customers
            for jj=1:length(midd)
                dd=dd+D(tour(ii+1),midd(jj));
            end
        end
        cluster_d(n)=dd/(number_of_customers*(length(midd)));
    end
end

cluster_d=max(cluster_d)-cluster_d+0.00001;
cluster_d(tour_index1)=0;
tour2_index=RouletteWheelSelection(cluster_d);
fragment2=s(index_of_depots(tour2_index):index_of_depots(tour2_index+1));
number_of_insertp=length(fragment2)-1;
[candidate,~]=unified_fast_feasibility_check(dem_delivery,dem_pickup,fragment1,fragment2,number_of_insertp,fragment_length,vehicle_cap,D);
if ~isempty(candidate) 
    orignal_Fitness=Fitness_evaluation([tour,fragment2(2:end)],1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
    nc=size(candidate,1);
    dom=zeros(1,nc);
    for nn=1:nc
        candidate_Fitness(nn,:)=Fitness_evaluation([xxx,candidate(nn,2:end)],1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
        dom(nn)=check_domination(candidate_Fitness(nn,:),orignal_Fitness);
    end
    candidate_Fitness(dom==-1,:)=[]; 
    candidate(dom==-1,:)=[]; 
    if ~isempty(candidate_Fitness) 
        [NC,~]=size(candidate_Fitness);
        [candidate,candidate_Fitness]=population_update(candidate,candidate_Fitness);
        selection_value=Nonlinear_modification(linspace(0,1,NC),gre); 
        selected_one=RouletteWheelSelection(selection_value);
        s(index_of_depots(tour_index1)+start_p-1:index_of_depots(tour_index1)+end_p-1)=inf;
        s=[s(1:index_of_depots(tour2_index)-1),candidate(selected_one,:),s(index_of_depots(tour2_index+1)+1:end)];
        s(s==inf)=[];
        individual=s;
        fitness=fitness-orignal_Fitness+candidate_Fitness(selected_one,:); 
        for a=1:length(individual)-1
            if individual(a)==1&&individual(a+1)==1
                individual(a)=inf;
            end
        end
        individual(individual==inf)=[];
        individual=[individual,zeros(1,N+vehicle_num-length(individual))];
    end
    if number_of_depots>sum(individual==1)
        fitness(1)=fitness(1)-2000;
    end
end
end

function [individual,fitness]=inter_TFR_exchange_SWAPstar(individual,fitness,dem_delivery,dem_pickup,customer_time,vehicle_cap,D,N,vehicle_num,gre,Miu,ALpha,BEta,max_fragment_length)
s=individual;
index_of_depots=find(s==1);
number_of_depots=length(index_of_depots)-1;
tour_index=randi(number_of_depots);
index_of_depots_t1=index_of_depots(tour_index);
tour=s(index_of_depots_t1:index_of_depots(tour_index+1)); 
number_of_customers=length(tour)-2;
cluster_d=zeros(1,number_of_depots);
for n=1:number_of_depots
    if n~=tour_index
        midd=s(index_of_depots(n)+1:index_of_depots(n+1)-1);
        dd=0;
        for ii=1:number_of_customers
            for jj=1:length(midd)-2
                dd=dd+D(tour(ii+1),midd(jj+1));
            end
        end
        cluster_d(n)=dd/(number_of_customers*(length(midd)-2));
    end
end
cluster_d=max(cluster_d)-cluster_d+0.00001;
cluster_d(tour_index)=0;
tour2_index=RouletteWheelSelection(cluster_d);
tours_index=[tour_index,tour2_index];
tours_index=sort(tours_index);
for i=1:2
    ti=tours_index(i);
    tours(i).s=s(index_of_depots(ti):index_of_depots(ti+1)); 
    tour=tours(i).s;
    number_of_customers=length(tour)-2;
    max_fragment_length=min(number_of_customers,max_fragment_length);
    fragment_length(i)=randi(max_fragment_length);
    start_p(i)=randi(number_of_customers-fragment_length(i)+1)+1;
    end_p(i)=start_p(i)+fragment_length(i)-1;
end
orignal_Fitness=Fitness_evaluation([tours(1).s,tours(2).s(2:end)],1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
fragment1=tours(1).s(start_p(1):end_p(1)); 
fragment2=tours(2).s; 
fragment2(start_p(2):end_p(2))=[];
number_of_insertp=length(fragment2)-1; 
[candidate1,~]=unified_fast_feasibility_check(dem_delivery,dem_pickup,fragment1,fragment2,number_of_insertp,fragment_length(1),vehicle_cap,D); 
fragment1=tours(2).s(start_p(2):end_p(2));
fragment2=tours(1).s;
fragment2(start_p(1):end_p(1))=[];
number_of_insertp=length(fragment2)-1; 
[candidate2,~]=unified_fast_feasibility_check(dem_delivery,dem_pickup,fragment1,fragment2,number_of_insertp,fragment_length(2),vehicle_cap,D);
if ~isempty(candidate1)&&~isempty(candidate2) 
    nc1=size(candidate1,1);
    candidate_fitness1=Fitness_evaluation(candidate1,nc1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
    [NC,~]=size(candidate_fitness1);
    [candidate1,candidate_fitness1]=population_update(candidate1,candidate_fitness1);
    selection_value=Nonlinear_modification(linspace(0,1,NC),gre);
    selected_one1=RouletteWheelSelection(selection_value);
    opt_candidate1=candidate1(selected_one1,:);

    nc2=size(candidate2,1);
    candidate_fitness2=Fitness_evaluation(candidate2,nc2,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
    [NC,~]=size(candidate_fitness2);
    [candidate2,candidate_fitness2]=population_update(candidate2,candidate_fitness2);
    selection_value=Nonlinear_modification(linspace(0,1,NC),gre);
    selected_one2=RouletteWheelSelection(selection_value);
    opt_candidate2=candidate2(selected_one2,:);

    s(index_of_depots(tours_index(1))+1:index_of_depots(tours_index(1)+1)-1)=inf;
    s(index_of_depots(tours_index(2))+1:index_of_depots(tours_index(2)+1)-1)=inf;
    s(s==0)=[];
    s=[s,opt_candidate1,opt_candidate2];
    individual=s;
    fitness=fitness-orignal_Fitness+candidate_fitness1(selected_one1,:)+candidate_fitness2(selected_one2,:);
    individual(individual==inf)=[];
    for a=1:length(individual)-1
        if individual(a)==1&&individual(a+1)==1
            individual(a)=inf;
        end
    end
    individual(individual==inf)=[];
    individual=[individual,zeros(1,N+vehicle_num-length(individual))];
end
end

function [individual,fitness]=intra_TFR(individual,fitness,dem_delivery,dem_pickup,customer_time,vehicle_cap,D,gre,Miu,ALpha,BEta,max_fragment_length)
s=individual;
s(s==0)=[];
index_of_depots=find(s==1);
number_of_depots=length(index_of_depots)-1;
tour_index=randi(number_of_depots);
tour=s(index_of_depots(tour_index):index_of_depots(tour_index+1)); 
number_of_customers=length(tour)-2;
max_fragment_length=min(number_of_customers,max_fragment_length);
fragment_length=randi(max_fragment_length); 
start_p=randi(number_of_customers-fragment_length+1)+1; 
end_p=start_p+fragment_length-1;
fragment1=tour(start_p:end_p);
fragment2=tour;
fragment2(start_p:end_p)=[]; 
number_of_insertp=length(fragment2)-1; 
[candidate,~]=unified_fast_feasibility_check(dem_delivery,dem_pickup,fragment1,fragment2,number_of_insertp,fragment_length,vehicle_cap,D);
if ~isempty(candidate) 
    orignal_Fitness=Fitness_evaluation(tour,1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta); 
    nc=size(candidate,1);
    dom=zeros(1,nc);
    candidate_Fitness=zeros(nc,2);
    for nn=1:nc
        candidate_Fitness(nn,:)=Fitness_evaluation(candidate(nn,:),1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
        dom(nn)=check_domination(candidate_Fitness(nn,:),orignal_Fitness);
    end
    candidate_Fitness(dom==-1,:)=[]; 
    candidate(dom==-1,:)=[];
    if ~isempty(candidate_Fitness) 
        [NC,~]=size(candidate_Fitness);
        [candidate,candidate_Fitness]=population_update(candidate,candidate_Fitness);
        selection_value=Nonlinear_modification(linspace(0,1,NC),gre);
        selected_one=RouletteWheelSelection(selection_value);
        individual(index_of_depots(tour_index):index_of_depots(tour_index+1))=candidate(selected_one,:);
        fitness=fitness-orignal_Fitness+candidate_Fitness(selected_one,:);
    end
end
end

function [individual,fitness]=intra_TFR_exchange(individual,fitness,dem_delivery,dem_pickup,customer_time,vehicle_cap,D,gre,Miu,ALpha,BEta,max_fragment_length)
s=individual;
index_of_depots=find(s==1);
number_of_depots=length(index_of_depots)-1;
tour_index=randi(number_of_depots);
tour=s(index_of_depots(tour_index):index_of_depots(tour_index+1));
orignal_Fitness=Fitness_evaluation(tour,1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta); 
number_of_customers=length(tour)-2;
while number_of_customers<2 
    tour_index=randi(number_of_depots);
    tour=s(index_of_depots(tour_index):index_of_depots(tour_index+1)); 
    number_of_customers=length(tour)-2;
end
max_fragment_length=min(number_of_customers-1,max_fragment_length);
for i=1:2
    fragment_length=randi(max_fragment_length); 
    start_p(i)=randi(number_of_customers-fragment_length+1)+1;
    end_p(i)=start_p(i)+fragment_length-1;
end
pp=sort([start_p,end_p]);
if pp(2)==pp(3)
    if rand>0.5
        if pp(2)-1>=2
            pp(2)=pp(2)-1;
            pp(1)=max(pp(1)-1,2);
        elseif pp(3)+1<=number_of_customers+1
            pp(3)=pp(3)+1;
            pp(4)=min(pp(4)+1,number_of_customers+1);
        end
    else
        if pp(3)+1<=number_of_customers+1
            pp(3)=pp(3)+1;
            pp(4)=min(pp(4)+1,number_of_customers+1);
        elseif pp(2)-1>=2
            pp(2)=pp(2)-1;
            pp(1)=max(pp(1)-1,2);
        end
    end
end
pp(1:2)=sort(pp(1:2));
pp(3:4)=sort(pp(3:4));
start_p(1)=pp(1);
end_p(1)=pp(2);
start_p(2)=pp(3);
end_p(2)=pp(4);
i=1;
candidate=[];
candidate_Fitness=[];
dom=[];
for nnn=1:4
    if nnn==2
        tour(start_p(1):end_p(1))=tour(end_p(1):-1:start_p(1));
    elseif nnn==3
        tour(start_p(2):end_p(2))=tour(end_p(2):-1:start_p(2));
    elseif nnn==4
        tour(start_p(1):end_p(1))=tour(end_p(1):-1:start_p(1));
    end
    candidate1_1=[tour(1:start_p(1)-1),tour(start_p(2):end_p(2)),tour(end_p(1)+1:start_p(2)-1),...
        tour(start_p(1):end_p(1)),tour(end_p(2)+1:end)];
    if Feasibility_check(candidate1_1,dem_delivery,dem_pickup,vehicle_cap)
        candidate(i,:)=candidate1_1;
        candidate_Fitness(i,:)=Fitness_evaluation(candidate(i,:),1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
        dom(i)=check_domination(candidate_Fitness(i,:),orignal_Fitness);
        i=i+1;
    end
end
if ~isempty(dom)
    candidate_Fitness(dom==-1,:)=[]; 
    candidate(dom==-1,:)=[];
end
if ~isempty(candidate_Fitness)
    [NC,~]=size(candidate_Fitness);
    [candidate,candidate_Fitness]=population_update(candidate,candidate_Fitness);
    selection_value=Nonlinear_modification(linspace(0,1,NC),gre); 
    selected_one=RouletteWheelSelection(selection_value);
    individual(index_of_depots(tour_index):index_of_depots(tour_index+1))=candidate(selected_one,:);
    fitness=fitness-orignal_Fitness+candidate_Fitness(selected_one,:);
end
end

function [Population,Fitness]=survival_of_the_fittest(Population,Fitness,popN)
[Population,Fitness]=population_update(Population,Fitness);
Population(popN+1:end,:)=[];
Fitness(popN+1:end,:)=[]; 
end

function [flag,rofinf,seq]=Feasibility_check(S,dem_delivery,dem_pickup,vehicle_cap)
flag=0;
index=find(S==1);
li=length(index);
seq=1;
for i=1:li-1
    s=S(index(i):index(i+1));
    remain_cap=vehicle_cap+sum(dem_delivery(s));
    remain_cap2=vehicle_cap-sum(dem_pickup(s));
    if remain_cap>=0&&remain_cap2>=0
        if sum((cumsum(dem_delivery(s)+dem_pickup(s))-remain_cap)<=0)==length(s)
            flag=1;
            rofinf=0;
        else
            flag=0;
            rofinf=2;
            break;
        end
    else
        flag=0;
        rofinf=1;
        break;
    end
end
end

function [candidate,rofinf]=unified_fast_feasibility_check(dem_delivery,dem_pickup,fragment1,fragment2,number_of_insertp,fragment_length,vehicle_cap,D)
rofinf=0; 
candidate=[];
L0=sum(dem_delivery([fragment1,fragment2]));
L1=sum(dem_pickup([fragment1,fragment2]));
if L0>vehicle_cap||L1>vehicle_cap
    rofinf=1;
else
    n_fload_fragment1=cumsum(dem_pickup(fragment1)-dem_delivery(fragment1));
    n_fload_fragment1_reverse=cumsum(dem_pickup(fragment1(end:-1:1))-dem_delivery(fragment1(end:-1:1)));
    J_fload_fragment1=n_fload_fragment1(end);
    J_fload_fragment1_reverse=n_fload_fragment1_reverse(end);
    J_lload_fragment1=sum(dem_delivery(fragment1));
    J_lload_fragment2=sum(dem_delivery(fragment2));
    n_fload_fragment2=[cumsum(dem_pickup(fragment2)-dem_delivery(fragment2))];
    n_fload_fragment2(end)=[];
    i=1;
    for ip=1:number_of_insertp
        ip_plus1=n_fload_fragment2(ip+1:end);
        if (sum(J_fload_fragment1<=vehicle_cap-L0-ip_plus1)==length(ip_plus1)&&...
                sum(J_lload_fragment1<=vehicle_cap-J_lload_fragment2-n_fload_fragment2(1:ip))==ip&&...
                sum(n_fload_fragment1<=vehicle_cap-L0-n_fload_fragment2(ip))==fragment_length)
                candidate(i,:)=[fragment2(1:ip),fragment1,fragment2(ip+1:end)];
                i=i+1;
        elseif (sum(J_fload_fragment1_reverse<=vehicle_cap-L0-ip_plus1)==length(ip_plus1)&&...
                sum(J_lload_fragment1<=vehicle_cap-J_lload_fragment2-n_fload_fragment2(1:ip))==ip&&...
                sum(n_fload_fragment1_reverse<=vehicle_cap-L0-n_fload_fragment2(ip))==fragment_length)
                candidate(i,:)=[fragment2(1:ip),fragment1(end:-1:1),fragment2(ip+1:end)];
                i=i+1;
        end
    end
end
end

function Unified_feasibility_check(S,customer_dem,vehicle_cap,customer_time,D,Miu,ALpha,BEta)
dem_delivery=customer_dem(:,1);
dem_pickup=customer_dem(:,2);
[Cap_feasibility,~,~]=Feasibility_check(S,dem_delivery,dem_pickup,vehicle_cap);
[~,~,Time_feasibility]=Fitness_evaluation(S,1,customer_time,D,dem_delivery,dem_pickup,Miu,ALpha,BEta);
if Cap_feasibility==1&&Time_feasibility==0

elseif Cap_feasibility==0&&Time_feasibility==0
    fprintf('Error:capacity-infeasible optimum!');
elseif Cap_feasibility==1&&Time_feasibility>0
    fprintf('Error:time-infeasible optimum!');
else
    fprintf('Error:capacity-time-infeasible optimum!');
end
end