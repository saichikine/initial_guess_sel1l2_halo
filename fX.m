function fX = fX(X, c, mu, Tmax)
	x = X(1);
 	y = X(2);
 	z = X(3);
 	xdot = X(4);
 	ydot = X(5);
 	zdot = X(6);
 	m = X(7);
 	ux = X(8);
 	uy = X(9);
 	uz = X(10);

	fX(1,1) = 0;
	fX(1,2) = 0;
	fX(1,3) = 0;
	fX(1,4) = 1;
	fX(1,5) = 0;
	fX(1,6) = 0;
	fX(1,7) = 0;
	fX(1,8) = 0;
	fX(1,9) = 0;
	fX(1,10) = 0;

	fX(2,1) = 0;
	fX(2,2) = 0;
	fX(2,3) = 0;
	fX(2,4) = 0;
	fX(2,5) = 1;
	fX(2,6) = 0;
	fX(2,7) = 0;
	fX(2,8) = 0;
	fX(2,9) = 0;
	fX(2,10) = 0;

	fX(3,1) = 0;
	fX(3,2) = 0;
	fX(3,3) = 0;
	fX(3,4) = 0;
	fX(3,5) = 0;
	fX(3,6) = 1;
	fX(3,7) = 0;
	fX(3,8) = 0;
	fX(3,9) = 0;
	fX(3,10) = 0;

	fX(4,1) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
	fX(4,2) = (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
	fX(4,3) = (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
	fX(4,4) = 0;
	fX(4,5) = 2;
	fX(4,6) = 0;
	fX(4,7) = -(Tmax*ux)/m^2;
	fX(4,8) = Tmax/m;
	fX(4,9) = 0;
	fX(4,10) = 0;

	fX(5,1) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
	fX(5,2) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
	fX(5,3) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
	fX(5,4) = -2;
	fX(5,5) = 0;
	fX(5,6) = 0;
	fX(5,7) = -(Tmax*uy)/m^2;
	fX(5,8) = 0;
	fX(5,9) = Tmax/m;
	fX(5,10) = 0;

	fX(6,1) = (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
	fX(6,2) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
	fX(6,3) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);
	fX(6,4) = 0;
	fX(6,5) = 0;
	fX(6,6) = 0;
	fX(6,7) = -(Tmax*uz)/m^2;
	fX(6,8) = 0;
	fX(6,9) = 0;
	fX(6,10) = Tmax/m;

	fX(7,1) = 0;
	fX(7,2) = 0;
	fX(7,3) = 0;
	fX(7,4) = 0;
	fX(7,5) = 0;
	fX(7,6) = 0;
	fX(7,7) = 0;
	fX(7,8) = -ux/(c*(ux^2 + uy^2 + uz^2 + 1/100000000)^(1/2));
	fX(7,9) = -uy/(c*(ux^2 + uy^2 + uz^2 + 1/100000000)^(1/2));
	fX(7,10) = -uz/(c*(ux^2 + uy^2 + uz^2 + 1/100000000)^(1/2));

	fX(8,1) = 0;
	fX(8,2) = 0;
	fX(8,3) = 0;
	fX(8,4) = 0;
	fX(8,5) = 0;
	fX(8,6) = 0;
	fX(8,7) = 0;
	fX(8,8) = 0;
	fX(8,9) = 0;
	fX(8,10) = 0;

	fX(9,1) = 0;
	fX(9,2) = 0;
	fX(9,3) = 0;
	fX(9,4) = 0;
	fX(9,5) = 0;
	fX(9,6) = 0;
	fX(9,7) = 0;
	fX(9,8) = 0;
	fX(9,9) = 0;
	fX(9,10) = 0;

	fX(10,1) = 0;
	fX(10,2) = 0;
	fX(10,3) = 0;
	fX(10,4) = 0;
	fX(10,5) = 0;
	fX(10,6) = 0;
	fX(10,7) = 0;
	fX(10,8) = 0;
	fX(10,9) = 0;
	fX(10,10) = 0;

end