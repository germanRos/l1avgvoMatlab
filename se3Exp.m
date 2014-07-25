function Rt = se3Exp(V)
	r = V(1:3);
	u = V(4:6);

	omega_x = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
	S = [omega_x  u'; 0 0 0 0];
	theta = norm(r);

	if theta == 0
		Rt = eye(4) + S;
	else
		Rt = eye(4) + S + ((1-cos(theta))/(theta*theta)).*(S*S) + ((theta-sin(theta))/(theta^3)).*(S*S*S);
	end
end
