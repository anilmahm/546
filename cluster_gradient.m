%calculate the sum of gradients in a cluster, the arrangement of clusters
%is passed as cluster, the measurements of each sensor as Y_l , L and c are
%defined in the paper, k is the iteration number , 
%function [cluster_index, grad_cluster] = cluster_gradient(x, k, cluster, num_cluster, cumulitive_cluster, Y_l, L, c) 
function [ grad_cluster] = cluster_gradient(x, k, cluster, num_cluster, cumulitive_cluster, Y_l, L, c) 
    
% find out which cluster to use by use of k 
    cluster_index = mod( k, num_cluster);
    if ( cluster_index == 0)
        cluster_index = num_cluster;
    end 
    
    %find the starting sensor 
    if( cluster_index == 1)  
        start = 1;
    else start = cumulitive_cluster(cluster_index-1) + 1; 
    end
    
    %sum up the gradients 
    
    grad_cluster = 0 ;
    for i = start: start+ cluster(cluster_index)-1
        grad_cluster = grad_cluster + function_prototype_grad( x, Y_l( i ), L , c) ;
    end 
           
end
