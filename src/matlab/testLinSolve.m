clear all
A = [1, 0 2; 0 1 2; 3 -3 2 ];
B = A-eye(3,3);
b = [1; 0;0];
restart = 100;
sol = gmres(B, b, restart, 1e-6, 100,[],[],[1/3; 1/3; 1/3]);

null(B)