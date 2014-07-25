function V = se3Log(Rt)
	epsilon = 10e-9;
	R = Rt(1:3, 1:3);
	T = Rt(1:3, 4);

	r = so3Log(R);
	nn = norm(r);
	Omega = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];

	if(nn < epsilon)
		t = T';
	else
		Vinv = eye(3) - 0.5*Omega + ((2.0*sin(nn) - nn - nn*cos(nn))/(2*nn*nn*sin(nn)))*Omega*Omega;
		t = (Vinv*T)';
	end

	V = [r, t];
end
