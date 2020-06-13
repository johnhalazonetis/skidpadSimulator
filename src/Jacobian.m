function [MatrixA, MatrixB] = Jacobian(v_x, v_y, delta, r, w, mu, kappa, Lf, Lr, m, Cdown, Cdrag, Wf, B, D, E, Theta, CarFriction, Iz)
    A11 = 0;
    A12 = v_x*cos(mu) - v_y*sin(mu);
    A13 = sin(mu);
    A14 = cos(mu);
    A15 = 0;
    
    A21 = (-kappa^2*(v_x*cos(mu) - v_y*sin(mu)))/(w*kappa - 1)^2;
    A22 = (-kappa*(v_y*cos(mu) + v_x*sin(mu)))/(kappa*w - 1);
    A23 = kappa*cos(mu)/(kappa*w - 1);
    A24 = -kappa*sin(mu)/(kappa*w - 1);
    A25 = 1;
    
    A31 = 0;
    A32 = 0;
    A33 = -(2*Cdrag*v_x + CarFriction/(v_x^2 + 1) - 2*sin(delta)*Cdown*D*Wf*v_x*sin((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D) + (sin(delta)*D^2*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*(E - 1)*(v_y + Lf)*r))/(v_x^2*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E*(v_y + Lf*r))/(v_x^2*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1)/(Theta + m);
    A34 = r + (sin(delta)*D^2*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*(E - 1))/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E)/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/((((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1)*(Theta + m));
    A35 = v_y + (sin(delta)*D^2*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*Lf*(E - 1))/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E*Lf)/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/((((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1)*(Theta + m));
    
    A41 = 0;
    A42 = 0;
    A43 = -r - (4*sin((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*Cdown*D*v_x*(Wf/2 - 1/2) + 2*cos(delta)*Cdown*D*Wf*v_x*sin((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D) + (2*cos((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*D^2*(Wf/2 - 1/2)*((981*m)/100 + Cdown*v_x^2)*((B*(E - 1)*(v_y - Lr*r))/(v_x^2*((v_y - Lr*r)^2/v_x^2 + 1)) - (B*E*(v_y - Lr*r))/(v_x^2*(B^2*(atan((v_y - Lr*r)/v_x))^2 + 1)*((v_y - Lr*r)^2/v_x^2 + 1))))/((E*(atan(B*atan((v_y - Lr*r)/v_x))) - B*(atan((v_y - Lr*r)/v_x))*(E - 1))^2 + 1) - (cos(delta)*D^2*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*(E - 1)*(v_y + Lf*r))/(v_x^2*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E*(v_y + Lf*r))/(v_x^2*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1))/m;
    A44 = ((2*cos((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*D^2*(Wf/2 - 1/2)*((B*(E - 1))/(v_x*((v_y - Lr*r)^2/v_x^2 + 1)) - (B*E)/(v_x*(B^2*(atan((v_y - Lr*r)/v_x))^2 + 1)*((v_y - Lr*r)^2/v_x^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/((E*(atan(B*atan((v_y - Lr*r)/v_x))) - B*(atan((v_y - Lr*r)/v_x))*(E - 1))^2 + 1) - (cos(delta)*D^2*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*(E - 1))/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E)/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1))/m;
    A45 = -v_x - ((2*cos((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*D^2*((B*Lr*(E - 1))/(v_x*((v_y - Lr*r)^2/v_x^2 + 1)) - (B*E*Lr)/(v_x*(B^2*(atan((v_y - Lr*r)/v_x))^2 + 1)*((v_y - Lr*r)^2/v_x^2 + 1)))*(Wf/2 - 1/2)*((981*m)/100 + Cdown*v_x^2))/((E*(atan(B*atan((v_y - Lr*r)/v_x))) - B*(atan((v_y - Lr*r)/v_x))*(E - 1))^2 + 1) + (cos(delta)*D^2*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*Lf*(E - 1))/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E*Lf)/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1))/m;
    
    A51 = 0;
    A52 = 0;
    A53 = ((4802*delta)/17 + 4*sin((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*Cdown*D*Lr*v_x*(Wf/2 - 1/2) - 2*cos(delta)*Cdown*D*Lf*Wf*v_x*sin((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D) + (2*cos((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*D^2*Lr*(Wf/2 - 1/2)*((981*m)/100 + Cdown*v_x^2)*((B*(E - 1)*(v_y - Lr*r))/(v_x^2*((v_y - Lr*r)^2/v_x^2 + 1)) - (B*E*(v_y - Lr*r))/(v_x^2*(B^2*(atan((v_y - Lr*r)/v_x))^2 + 1)*((v_y - Lr*r)^2/v_x^2 + 1))))/((E*(atan(B*atan((v_y - Lr*r)/v_x))) - B*(atan((v_y - Lr*r)/v_x))*(E - 1))^2 + 1) + (cos(delta)*D^2*Lf*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*(E - 1)*(v_y + Lf*r))/(v_x^2*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E*(v_y + Lf*r))/(v_x^2*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1))/Iz;
    A54 = -((2*cos((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*D^2*Lr*(Wf/2 - 1/2)*((B*(E - 1))/(v_x*((v_y - Lr*r)^2/v_x^2 + 1)) - (B*E)/(v_x*(B^2*(atan((v_y - Lr*r)/v_x))^2 + 1)*((v_y - Lr*r)^2/v_x^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/((E*(atan(B*atan((v_y - Lr*r)/v_x))) - B*(atan((v_y - Lr*r)/v_x))*(E - 1))^2 + 1) + (cos(delta)*D^2*Lf*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*(E - 1))/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E)/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1))/Iz;
    A55 = -((cos(delta)*D^2*Lf*Wf*cos((atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*D)*((B*Lf*(E - 1))/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)) - (B*E*Lf)/(v_x*((v_y + Lf*r)^2/v_x^2 + 1)*(B^2*(delta - (atan((v_y + Lf*r)/v_x)))^2 + 1)))*((981*m)/100 + Cdown*v_x^2))/(((atan(B*(delta - atan((v_y + Lf*r)/v_x))))*E - B*(E - 1)*(delta - (atan((v_y + Lf*r)/v_x))))^2 + 1) - (2*cos((atan(E*atan(B*atan((v_y - Lr*r)/v_x)) - B*atan((v_y - Lr*r)/v_x)*(E - 1)))*D)*D^2*Lr*((B*Lr*(E - 1))/(v_x*((v_y - Lr*r)^2/v_x^2 + 1)) - (B*E*Lr)/(v_x*(B^2*(atan((v_y - Lr*r)/v_x))^2 + 1)*((v_y - Lr*r)^2/v_x^2 + 1)))*(Wf/2 - 1/2)*((981*m)/100 + Cdown*v_x^2))/((E*(atan(B*atan((v_y - Lr*r)/v_x))) - B*(atan((v_y - Lr*r)/v_x))*(E - 1))^2 + 1) + 21109/50)/Iz;
    
    MatrixA = [A11 A12 A13 A14 A15 ; A21 A22 A23 A24 A25 ; A31 A32 A33 A34 A35 ; A41 A42 A43 A44 A45 ; A51 A52 A53 A54 A55];
    
    B11 = 0;
    B12 = 0;
    
    B21 = 0;
    B22 = 0;
    
    B31 = 2*20.2086*m/(Theta + m);
    B32 = (D*Wf*sin(D*atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*cos(delta)*(Cdown*v_x^2 + (981*m)/100) - (D^2*Wf*cos(D*atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*sin(delta)*(B*(E - 1) - (B*E)/(B^2*(delta - atan((v_y + Lf*r)/v_x))^2 + 1))*(Cdown*v_x^2 + (981*m)/100))/((E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x)))^2 + 1))/(Theta + m);
    
    B41 = 0;
    B42 = (D*Wf*sin(D*atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*sin(delta)*(Cdown*v_x^2 + (981*m)/100) + (D^2*Wf*cos(D*atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*cos(delta)*(B*(E - 1) - (B*E)/(B^2*(delta - atan((v_y + Lf*r)/v_x))^2 + 1))*(Cdown*v_x^2 + (981*m)/100))/((E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x)))^2 + 1))/m;
    
    B51 = 0;
    B52 = ((4802*v_x)/17 + D*Lf*Wf*sin(D*atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*sin(delta)*(Cdown*v_x^2 + (981*m)/100) + (D^2*Lf*Wf*cos(D*atan(E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x))))*cos(delta)*(B*(E - 1) - (B*E)/(B^2*(delta - atan((v_y + Lf*r)/v_x))^2 + 1))*(Cdown*v_x^2 + (981*m)/100))/((E*atan(B*(delta - atan((v_y + Lf*r)/v_x))) - B*(E - 1)*(delta - atan((v_y + Lf*r)/v_x)))^2 + 1))/Iz;
    
    MatrixB = [B11 B12; B21 B22; B31 B32; B41 B42; B51 B52];
    
    
end





