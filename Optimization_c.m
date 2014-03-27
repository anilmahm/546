%
%   Implementation IAG, IG based on Hero Blatt Gauchman Paper 
%
%

close all 
clear all 

L = 50 ; % number of sensors 
L_orig = L ; % for using in equation 1 of ICAG 
c = 10 ; % variable c used in paper 
mean1 = 10; % mean used in paper ( both were 10 ) 

sigma1 = 1;  % variance used in paper 
sigma2 = sqrt(10);% variance used in paper 
%generate sensor measurements , first half with mu = 10, sigma = 1
first_half = normrnd( mean1, sigma1, [1,floor(L/2)]);
second_half = normrnd( mean1, sigma2, [1,L - floor(L/2)]);

% set to 1 to cause damage to sensors, see code for ICAG 
destroy_sensors  = 0;

%create the measurements according to paper Hero Blatt Gauchman , 50%
%measurements had high variance, 50% had low variance
Y_l = [first_half  second_half];

% save('Y_l.mat', 'Y_l')
load('Y_l.mat','Y_l');
%Y_l(1)=0;


%max number of iterations 
MAX_ITER = 1000;

% the iterates, for IAG of Hero blatt gauchman , SIG is the standard
% incremnetal gradient method, DAG is the method I implemented named ICAG
% in paper 
x_val_IAG = zeros(1, MAX_ITER);
x_val_SIG = zeros(1, MAX_ITER);
x_val_DAG = zeros(1, MAX_ITER);

%backup of gradients, used to form the aggregate gradient 
temp_grad_IAG = zeros(1, MAX_ITER + L );
temp_grad_SIG = zeros(1, MAX_ITER + L );
temp_grad_DAG = zeros(1, MAX_ITER + L );
%step size

%optimal step sizes found , only the last one is effective, commented out
%when wanted to use the previous one these stepsizes are optimal for IAG
%and ICAG respectively for the function Fair used in paper 
Mu =0.95; %447
Mu =4.78; %  89
%Mu = 5

%these are the optimal step sizes for sum of quadratics problem 
% Mu =0.37;  %371
% Mu = 1.86;  %73


%wanted to find the optimal point correct upto 4 decimal places 
epsilon = 0.00005;

%initial point 
x_init = 0;



%%%%
%%%%            Incrimental clustered AGGREGATE GRADIENT METHOD   (ICAG) 
%%%%



%cluster = [ 5  15 10  8  12 ];
%cluster = [ 15 15 ];

% 50 sensors, 10 clusters of size 5 each 
cluster = [ 5  5 5 5 5 5 5 5 5 5 ];

% for test purposes 
%cluster = [ 1  9 1 9 1 9 1 9 1 9 ];
% cluster = [ 1  9 1 9 1 9 1 9 1 9 ];
% cluster = [ 3  7 3 7 3 7 3 7 3 7 ];

%number of clusters 
num_cluster = length(cluster);

% cumulitive sum of sensors in cluster , needed for ICAG method that i
% implemented 
cumulitive_cluster = cumsum(cluster); 


%initialization phase
x_val_DAG(1) = x_init ;
temp_grad_DAG(1) = cluster_gradient(x_val_DAG(1),1, cluster, num_cluster, cumulitive_cluster, Y_l, L , c);
d= temp_grad_DAG(1);

%initialization of N_c points as mentioned in paper 
for k = 1:num_cluster-1
    x_val_DAG(k+1) = x_val_DAG(k) - Mu * (1/cumulitive_cluster(k))* d;
   % x_val_DAG(k+1) = x_init;
   temp_grad_DAG(k+1) = cluster_gradient(x_val_DAG(k+1),k+1, cluster, num_cluster, cumulitive_cluster, Y_l, L , c);
    d = d + temp_grad_DAG(k+1);
end

start = k+1;
%optimization iterations
grad_count = start;


for k = start: MAX_ITER   
    

% uncomment to permute the clusters , for reclustering at iteration 30 
%     if(k ==30)
%         Y_l = Y_l( randperm(length(Y_l)));
%     end
    
    if( destroy_sensors ==1)
        %DAMAGE TO SENSORS
        
        if ( k == 30 )
            
            
            damaged_cluster = [];
           % damaged_cluster = [1];   % fill this up with cluster number to totally damage cluster , not supported without central processing system but implemented for testing purpose
            % damage the sensors, of clusters of size 5, maintain
            % corresponxdence with damaged cluster, the commented line
            % damages one cluster 
            cluster = [ 3 4 5 2 5 5 4 1 4 2 ];
           %cluster = [ 4 5 2 5 5 4 1 4 2 ];  % (first cluster is damaged, if
           %damaged_cluster = [1];
            % remove the sensor measurements based on sensor damage this
            % has to be done depending on the damage caused by the previous
            % two lines manually, 
           Y_l = [Y_l(1:3)  Y_l(6:9)  Y_l(11:15)  Y_l(16:17)  Y_l(21:30) Y_l(31:34) Y_l(36) Y_l(41:44) Y_l(46:47) ];
           
           % use this if cluster 1 is damaged , 
       %   Y_l = [  Y_l(6:9)  Y_l(11:15)  Y_l(16:17)  Y_l(21:30) Y_l(31:34) Y_l(36) Y_l(41:44) Y_l(46:47) ];
          
            
            cluster_index = mod( grad_count, num_cluster);
            if ( cluster_index == 0)
                cluster_index = num_cluster;
            end
            
            damaged_cluster = sort(damaged_cluster);
            
            
            % this loop is only needed for handling the damage of a whole
            % cluster, is not related to any of the algorithms given in
            % paper, just implemented to teest how things would go , to
            % test one has to fill the array damaged_cluster with some
            % cluster ids ( not important to understand or test paper that
            % i submitted)
            
            %find out which parts of the aggregate gradient has to be
            %removed, and update cluster number 
            for cluster_damage_count = 1: length(damaged_cluster)
                
               
                if( damaged_cluster(cluster_damage_count) < cluster_index )
                    %damaged_shift = cluster_index - damaged_cluster;
                    damaged_shift = cluster_index - damaged_cluster(cluster_damage_count) ;
                else
                    % damaged_shift = cluster_index - 1 + num_cluster - damaged_cluster;
                    damaged_shift = cluster_index - 1 + num_cluster - damaged_cluster(cluster_damage_count);
                end
                
                %could use k or grad_count
                d = d - temp_grad_DAG( grad_count - damaged_shift) ;
                
                temp_grad_DAG(grad_count-damaged_shift: grad_count) = temp_grad_DAG(grad_count-damaged_shift +1: grad_count+1);
                grad_count = grad_count -1 ;
                num_cluster = num_cluster-1;
            end
            
            
            %update number of sensors and cumulitive sum of sensor numbers 
            
            L = sum(cluster);
            
            cumulitive_cluster = cumsum(cluster);
            
            
        end
        
    end
    
    
    %equation 1 of ICAG to update and find iterate 
    x_val_DAG(k+1) = x_val_DAG(k) - Mu* (1/L_orig) * d;
    
    %calculate the sum of gradients for this cluster 
    temp_grad_DAG(k+1) = cluster_gradient(x_val_DAG(k+1),k+1, cluster, num_cluster, cumulitive_cluster, Y_l, L , c);
    %  d = d - temp_grad_DAG(k+1-num_cluster) + temp_grad_DAG(k+1);
    
    %update direction by using equation 2 of ICAG remove previous sum of gradient
    %of this cluster and add new sum of gradients for this cluster
    d = d - temp_grad_DAG(grad_count+1-num_cluster) + temp_grad_DAG(grad_count+1);
    %we need a separate counter for gradient to handle sensor damage. 
    grad_count = grad_count+1;
end






%%%%%
%%%%%            INCREMENTAL AGGREGATE GRADIENT METHOD
%%%%%

%initialization phase 
x_val_IAG(1) = x_init ;
temp_grad_IAG(1) = function_prototype_grad( x_val_IAG(1), Y_l(1), L, c) ;
d= temp_grad_IAG(1);

for k = 1:L-1
  x_val_IAG(k+1) = x_val_IAG(k) - Mu* (1/k) * d;
  %x_val_IAG(k+1)= x_init;
 func_index =  mod( k+1, L );
    if ( func_index == 0)
        func_index = L;
    end
    temp_grad_IAG(k+1) = function_prototype_grad( x_val_IAG(k+1), Y_l( func_index ), L , c) ;
    d = d + temp_grad_IAG(k+1);
end

start = k+1;
%optimization iterations

for k = start: MAX_ITER
    x_val_IAG(k+1) = x_val_IAG(k) - Mu* (1/L) * d;
    func_index =  mod( k+1, L );
    if ( func_index == 0)
        func_index = L;
    end
    temp_grad_IAG(k+1) = function_prototype_grad( x_val_IAG(k+1), Y_l( func_index ), L , c) ;
    d = d - temp_grad_IAG(k+1-L) + temp_grad_IAG(k+1);
end





%%%%
%%%%         STANDARD INCREMENTAL GRADIENT METHOD
%%%%

%initialization phase and optimization phase are in the same loop , as it
%does not requrie L or N_c initial points 

x_val_SIG(1) = x_init ;
temp_grad_SIG(1)= function_prototype_grad( x_val_SIG(1), Y_l(1), L, c) ;
d = temp_grad_SIG(1);


for k = 1:MAX_ITER
    x_val_SIG(k+1) = x_val_SIG(k) - Mu*  d;
    func_index =  mod( k+1, L );
    if ( func_index == 0)
        func_index = L;
    end
    temp_grad_SIG(k+1) = function_prototype_grad( x_val_SIG(k+1), Y_l( func_index ), L , c) ;
    d=  temp_grad_SIG(k+1);
end




% find the index of the point where the algorithm converged to 4 decimal
% places, , here I was checking convergence of the ICAG algorithm, modified
% this code when i needed to check for IAG , IG or standard incremental one
% never converges as was shown in paper , it hits a limit cycle 

res = x_val_DAG(MAX_ITER);

for conv_count = MAX_ITER : -1: 1
    if (abs(res - x_val_DAG(conv_count) ) > epsilon)
        break ;
    end 
end 

conv_count
x_val_DAG(1001)
    










figure
plot(1:MAX_ITER+1,x_val_SIG,1:MAX_ITER+1,x_val_IAG,1:MAX_ITER+1, x_val_DAG, 'LineWidth',1.4)
hleg1 = legend('IG','IAG', 'ICAG',...
                'Location','SouthEast' );
print -depsc plot.eps


