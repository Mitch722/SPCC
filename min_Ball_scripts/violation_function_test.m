% test violation function
clear all

R = 1;
c = [0,0];
global_lite = [1 0; 0 1; 1 1; 2 2 ];

[ no_violate, points, viol_fact ] = violation_function( R, c, global_lite );

