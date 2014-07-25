function R = so3Exp(V)
	theta = norm(V);
	if(theta == 0)
		R = eye(3);
	else
		v = V ./ theta;
		Omega = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
		R = eye(3) + sin(theta).*Omega + (1-cos(theta)).*(Omega*Omega);
	end
end
