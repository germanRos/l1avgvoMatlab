function J = derivativeCardan(v)
	ca1 = cos(v(1));
	sa1 = sin(v(1));
	ca2 = cos(v(2));
	sa2 = sin(v(2));
	ca3 = cos(v(3));
	sa3 = sin(v(3));	

	J1 = [	0,			0,				0,		0;
		ca1*ca3*sa2 - sa1*sa3,	-ca3*sa1 - ca1*sa2*sa3,		-ca1*ca2,	0;
		ca3*sa1*sa2 + ca1*sa3,	ca1*ca3 - sa1*sa2*sa3,		-ca2*sa1,	0;
		0,			0,				0,		0
	     ];

	J2 = [	-ca3*sa2,		sa2*sa3,			ca2,		0;
		ca2*ca3*sa1,		-ca2*sa1*sa3,			sa1*sa2,	0;
		-ca1*ca2*ca3,		ca1*ca2*sa3,			-ca1*sa2,	0;
		0,			0,				0,		0
	     ];

	J3 = [	-ca2*sa3,		-ca2*ca3,			0,		0;
		ca1*ca3 - sa1*sa2*sa3,	-ca3*sa1*sa2 - ca1*sa3,		0,		0;
		ca3*sa1 + ca1*sa2*sa3,	ca1*ca3*sa2 - sa1*sa3,		0,		0;
		0,			0,				0,		0
	     ];
	
	J = {J1, J2, J3};
end
