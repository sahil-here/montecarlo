% Takes dimension dim, Number of Runs R and Number of sample points N as input
% Returns the estimated volume in terms of MEAN, Standard deviation STD, and error bar ERR
% Everthing is returned in 4 digit precision
function [MEAN, STD, ERR] = montecarlo(dim, R, N)
    % the volume of a cube of edge length 2
    cubevol = 2^dim;
    S = zeros(R,1);
    for i = 1:R
        % draw the samples in the range of (-1,1)
        X = 2*rand(N,dim) - 1;
        % count the number of points inside the given volume
        C1 = sum(isInside(X));
        % compute the estimated volume
        V = ((C1)/(N))*cubevol;
        S(i) = V;
    end
    MEAN = round(mean(S),3);
    STD = round(std(S),3);
    ERR = round((1.96*std(S)/sqrt(R)),3);
end

%Takes a vector X and check which all points in the vector are inside the hypersphere
function i = isInside(X)
    R = sqrt(sum(X.^2,2));
    i = (R <= 1);
end

% Takes dimension dim, Number of Runs R and Number of sample points N as input
% Returns the estimated volume in terms of MEAN, Standard deviation STD, and error bar ERR
% Everthing is returned in 4 digit precision
% Variance reduction technique - Antithetic Sampling
function [MEAN, STD, ERR] = montecarloWithVarRed(dim, R, N)
    cubevol = 2^dim;
    S = zeros(R,1);
    for i = 1:R
        % draw the samples in the range of (-1,1)
        X = 2*rand(N,dim) - 1;
        Y = X*(-1);
        % count the number of points inside the given volume
        C1 = sum(isInside(X));
        C2 = sum(isInside(Y));
        % compute the estimated volume
        V = ((C1+C2)/(2*N))*cubevol;
        S(i) = V;
    end
    MEAN = round(mean(S),3);
    STD = round(std(S),3);
    ERR = round((1.96*std(S)/sqrt(R)),3);
end

% Takes dimension dim, Number of Runs R and Number of sample points N as input
% Returns the estimated volume in terms of MEAN, Standard deviation STD, and error bar ERR
% Everthing is returned in 4 digit precision
% Variance reduction technique - Randomized Quasi Monte Carlo Method
function [MEAN, STD, ERR] = randQuasiMC(dim, R, N)
    % the volume of a cube of edge length 2
    cubevol = 2^dim;
    N1 = N/20;
    S = zeros(R,1);
    for i = 1:R
        points = 0;
        %Here I am splitting to avoid memory issues for higher dimensions
        for n = 1:20
            p = sobolset(dim);
            % draw the random numbers
            Y = rand(N1,dim);
            % draw the numbers from sobol set
            X1 = net(p,N1);
            % create randomized quasi numbers
            X = 2*mod((X1+Y),1)-1;
            % count the number of points inside the volume
            CX1 = sum(isInside(X));
            points = points + CX1;
        end
        V = ((points)/(20*N1))*cubevol;
        S(i) = V;
    end
    MEAN = round(mean(S),3);
    STD = round(std(S),3);
    ERR = round((1.96*std(S)/sqrt(R)),3);
end

% function which returns the estimated volume by Cube or Grid Based Integration
%Takes K - number of intervals and dimension dim as inputs
function vol = cubebased(K, dim)
    length = (2/K);
    inside = 0;
    outside = 0;
    bisect = 0;
    MAXROWS=dim;
    MAXVALUES=K;
    points=linspace(-1,1,K+1);
    arrs = ones(1,MAXROWS);
    status = false;
    svol = length^dim;
    while(status==false)
        %test for exit condition i.e. we should stop when we get the last point
        %K,K,K
        if(sum(arrs) == MAXVALUES*MAXROWS)
            status = true;
        end
        %we get the coordinates here
        coordinates = zeros(1,MAXROWS);
        for j=1:MAXROWS
            %construct coordinates of the cube's center point
            coordinates(j) = points(arrs(j));
        end
        %check whether this cube is inside/outside/bisected here
        coordinates = coordinates + length/2; %to get the centroid
        boxtype = findBoxType(coordinates, K, dim);
        if(boxtype == 1) %1 means inside
            inside = inside + 1;
        elseif(boxtype == 2) %2 means outside
            outside = outside +1;
        else
            bisect = bisect + 1;
        end

        %increment loop variable
        change = true;
        %starting from innermost loop
        i = MAXROWS;
        while(change==true && i>=1)
            %increment the innermost variable and check if spill overs
            arrs(i) = arrs(i) + 1;
            if(arrs(i)>MAXVALUES)
                arrs(i) = 1; %reinitialize loop variable
                %change the upper variable by one
                change = true;
            else
                change = false; %stop as the upper levels of the loop are unaffected
            end
            %move to upper level of the loop
            i = i - 1;  
        end
    end
    vol = (inside+(bisect/2))*svol;
end

%function which returns 1 if the cube is inside, 2 for outside and 3 if it intersects the hypersphere 
function val = findBoxType(coordinates, K, dim)
    length = (2/K);
    radius = 1;
    diagonal = (sqrt(dim)*length)/2;
    distance = sqrt(sum(coordinates.^2,2));
    if(distance<=(radius-diagonal))
        val = 1; %for inside
    elseif(distance>=(radius+diagonal))
        val = 2; %for outside
    else
        val = 3; %for intersect
    end
end