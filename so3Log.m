function V = so3Log(R)
	theta = acos(0.5*(trace(R)-1));
	omega = (1/(2*sin(theta))) .* [R(3, 2) - R(2, 3), R(1, 3) - R(3, 1), R(2, 1) - R(1, 2)];

	if(theta == 0)
		V = zeros(1, 3);
	else
		V = theta.*omega;
	end
end
