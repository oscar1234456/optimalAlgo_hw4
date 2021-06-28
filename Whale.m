%Author: Oscar Chen
%Title: The Whale Algorithm
searchAgentsNum = 50; %The Max number of the whales
Max_Iteration = 50; %The Max number of iterations
b = 1; %The constant of exp function (Spiral)

dim = 2; %The dimension of the problem

ub = 600; %The upper(&lower) bound of the problem
lb = -600;

FitnessBest = 99999999; %record the best fitness

result = [];

%------Initialize the Whales-------
whalesPositions = rand(searchAgentsNum, dim).*(ub-lb)+lb; 


for i=1:size(whalesPositions,1)
    count = 0;
    x1 = whalesPositions(i,1);
    x2 = whalesPositions(i,2);
    fitness = (((x1^2) + (x2^2))/4000) - (cos(x1)*cos((x2)/sqrt(2)))+1; %Fitness Function
    if fitness < FitnessBest
        FitnessBest = fitness;  %update now best fitness value
        BestSearchAgent = whalesPositions(i,:);
    end
end

for t=1:Max_Iteration
    a = 2-(t-1)*((2)/Max_Iteration); %Caculate a (linearly descend)
    for i=1:size(whalesPositions,1)
        r1 = rand();
        r2 = rand();
        
        A = 2*a*r1-a;
        C = 2*r2;
       
        l = (a)*rand-1; %The coefficient of exp function
        p = rand();
        
        if p<0.5
            if abs(A)>=1
                %--------Search for prey(exploration phase)--------
                randomSearchAgentIndex = floor(searchAgentsNum*rand()+1);
                randomSearchAgent = whalesPositions(randomSearchAgentIndex,:);
                D_randomSearchAgent = abs((C*randomSearchAgent) - whalesPositions(i,:));
                whalesPositions(i,:) = randomSearchAgent - (A*D_randomSearchAgent);
            elseif abs(A)<1
                %--------Shrinking encircling mechanism---------
                D_BestSearchAgent = abs(C*BestSearchAgent - whalesPositions(i,:));
                whalesPositions(i,:) = BestSearchAgent - A*D_BestSearchAgent;
            end
        elseif p>=0.5
            %-------Spiral updating position--------
            distance2Leader = abs(BestSearchAgent-whalesPositions(i,:));
            whalesPositions(i,:) = distance2Leader*exp(b*l)*cos(2*pi*l)+BestSearchAgent;
        end    
    end
    for m=1:size(whalesPositions,1)
        x1 = whalesPositions(m,1);
        x2 = whalesPositions(m,2);
        fitness = (((x1^2) + (x2^2))/4000) - (cos(x1)*cos((x2)/sqrt(2)))+1; %Fitness Function
        if fitness < FitnessBest
            FitnessBest = fitness;
            BestSearchAgent = whalesPositions(m,:);%update now best fitness value
        end
    end
    result(t) = FitnessBest;
end
nexttile
title('每代目標函數最小值')
plot(1:Max_Iteration, result);
xlabel('Generation') 
ylabel('min Object Function Value') 

