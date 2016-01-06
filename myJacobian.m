function [q] = myJacobian(dh)
%% dh matrix should have parameters in order alpha, ai, di and theta 
%% for theta put 1 if joint is revolute else 0 if joint is prismatic
thetas = sym('theta', [1, size(dh,1)]);
%R = sym('r',[size(dh,1),(3 ,3)]);
J = sym('j',[6,size(dh,1)]);
T = eye(4);
% for i = 1 : (size(dh,1)+1)
%     R(1:3,1:3,i) = [0 0 1; 0 1 0; 0 0 1];
% end
R = sym('r',[3,3]);
R = eye(3);
digits(2);
for i = 1:size(dh,1)
    alpha = dh(i,1);
    ai = dh(i,2);
    di = dh(i,3);
    A(1,1,i)= cos(thetas(1,i));
    A(1,2,i)= -sin(thetas(1,i))*cos(alpha);
    A(1,3,i)= sin(thetas(1,i))*sin(alpha);
    A(1,4,i)= ai*cos(thetas(1,i));
    A(2,1,i)= sin(thetas(1,i));
    A(2,2,i)= cos(thetas(1,i))*cos(alpha)
    A(2,3,i)= -cos(thetas(1,i))*sin(alpha);
    A(2,4,i)= ai*sin(thetas(1,i));
    A(3,1,i)= 0;
    A(3,2,i)= sin(alpha);
    A(3,3,i)= cos(alpha);
    A(3,4,i)= di;
    A(4,1,i)= 0;
    A(4,2,i)= 0;
    A(4,3,i)= 0;
    A(4,4,i)= 1;
    T = T* A(:,:,i);
    digits(2);
    R= vpa(R*A(1:3,1:3,i));
    k = R(1:3,1:3)*[0; 0 ; 1];
    if dh(i,4)==1
       J(4,i) = k(1,1);
       J(5,i) = k(2,1);
       J(6,i) = k(3,1);
    else
        J(4:6,i)=0;
    end
    
end
digits(2);
for i = 1 : size(dh,1)
    J(1,i) = diff(vpa(T(1,4)),thetas(1,i));
    J(2,i) = diff(vpa(T(2,4)),thetas(1,i));
    J(3,i) = diff(vpa(T(3,4)),thetas(1,i)); 

end    
q = J;
disp(J);
end
