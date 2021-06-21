SearchAgents_no = 30;
Max_iteration = 200;
dim = 2;
ub = 600;
lb = -600;

leader_score = 99999999;

Positions = rand(SearchAgents_no, dim).*(ub-lb)+lb;

result = [];

for i=1:size(Positions,1)
    count = 0;
%     for j=1:size(Positions,2)
%         %count = count+(Positions(i,j)^2); 
%     end
    x1 = Positions(i,1);
    x2 = Positions(i,2);
    fitness = (((x1^2) + (x2^2))/4000) - (cos(x1)*cos((x2)/sqrt(2)))+1;
    if fitness < leader_score
        leader_score = fitness;
        leader_pos = Positions(i,:);
    end
end

for t=1:Max_iteration
    a = 2-(t-1)*((2)/Max_iteration);
    a2 = -1+(t-1)*((-1)/Max_iteration);
    for i=1:size(Positions,1)
        r1 = rand();
        r2 = rand();
        
        A = 2*a*r1-a;
        C = 2*r2;
        
        b = 1;
        l = (a2-1)*rand+1;
        p = rand();
        
        if p<0.5
            if abs(A)>=1
                rand_leader_index = floor(SearchAgents_no*rand()+1);
                X_rand = Positions(rand_leader_index,:);
                D_X_rand = abs((C*X_rand) - Positions(i,:));
                Positions(i,:) = X_rand - (A*D_X_rand);
            elseif abs(A)<1
                D_Leader = abs(C*leader_pos - Positions(i,:));
                Positions(i,:) = leader_pos - A*D_Leader;
            end
        elseif p>=0.5
            distance2Leader = abs(leader_pos-Positions(i,:));
            Positions(i,:) = distance2Leader*exp(b*l)*cos(2*pi*l)+leader_pos;
        end    
    end
    for m=1:size(Positions,1)
%     count = 0;
%         for y=1:size(Positions,2)
%             count = count+(Positions(m,y)^2); 
%         end
        x1 = Positions(m,1);
        x2 = Positions(m,2);
        fitness = (((x1^2) + (x2^2))/4000) - (cos(x1)*cos((x2)/sqrt(2)))+1;
        if fitness < leader_score
            leader_score = fitness;
            leader_pos = Positions(m,:);
        end
    end
    result(t) = leader_score;
end
nexttile
title('每代染色體Fitness最大值')
plot(1:Max_iteration, result);
xlabel('Generation') 
ylabel('min Fitness') 

