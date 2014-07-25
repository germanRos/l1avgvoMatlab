% L1AVGVO v1.0 by German Ros (gros@cvc.uab.es)
% This code performs Robust Stereo Visual Odometry from a set of
% matches between two views provided in Data and a calibration structure "calib".
%
% 	Input:
%		Data: 		Nx8 structure containing point matches (u, v) between two stereo views
%				i.e., (left, right)_tk+1 (left, right)_tk
%
%		calib:  	A structure containing calib.K a 3x3 matrix of intrinsic parameters
%				and calib.B a scalar that specifies the baseline in meters
%
%		options:	A structure with all the parameters required to run the method:
%			* [options.nModels] Number of putative models generated from the input data (e.g., 100 - 1000)
%			* [options.maxIters] Maximum amount of iterations for the iterative optimization of L1 averaging (e.g., 100- 1000)
%			* [options.threshold] Projection threshold for an optional outlier refinement step (e.g., 0.5 - 3.5)	
%			* [options.refinement] In case you want to activate the optional outlier refinement (true, false)
%			* [options.subModels] Number of models selected for the averaging (the top K models) (e.g., 500)
%			* [options.epsilon] Used in the stop condition of the L1 averaging (e.g., 1e-6)
%
%
%	Output:
%		sol:	A relative pose [R, T] in SE(3) as a 3x4 matrix
%
%
% This is a Matlab version of a C++ code available at https://github.com/germanRos/l1avgvo
% If you find any problem I would really appreciate if you inform me by e-mail.
% Of course, the C++ implementation is much faster, so in case you want to compare times, please
% always refer to the C++ version.
% 
% If you plan to use this code in a paper, please remember to cite the following paper:
%	@inproceedings{ros13,
% 	author = {Ros, G. and Guerrero, J. and Sappa, A. D. and Ponsa, D. and L\'{o}pez, A. M.},
% 	title = {Fast and Robust l1-averaging-based Pose Estimation for Driving Scenarios},
% 	booktitle = {Proceedings of the British Machine Vision Conference},
% 	year = {2013},
% 	address = {Bristol, UK},
%	}
%
%
function [sol] = L1AVGVO(Data, calib, options, ~)
	% generate N models and the Reduced Meassurement Matrices (l2-embbedings)
	[models, RMM_left, RMM_right] = generateModels(Data, calib.K, calib.B, options);
	
	% we estimate their scores
	[scores] = estimateScores(models, RMM_left, RMM_right);
	
	% select the N-"best" according to our metric
	NM = min(options.subModels, size(models, 2));   
   	indicesBestModels = selectModels(scores, NM);
    	bestModels = models(indicesBestModels);
    
	% model averaging
	sol = l1averaging(bestModels, options.epsilon, options.maxIters);

	%optional refienment
	if(options.refinement)
		sol = recalculateModelInliers(sol, X_left_p, pts_left_c, pts_right_c, eye(3), B, options.threshold);
	end
end

function S = l2averaging(Models, epsilon, maxIters)
	Tr = Models{1};
	Tr = [Tr; 0 0 0 1];

	N = numel(Models);

	done = false;
	iters = 1;
	while(~done && (iters < maxIters) )
		r = zeros(1, 6);
		for i=1:N
			Mi = Models{i};
			Tri = [Mi; 0 0 0 1];
			vi = se3Log(inv(Tr)*Tri);
			r = r + vi;
		end
		r = r ./ N;

		if(norm(r) < epsilon)
			done = true;
		else
			Tr = Tr * se3Exp(r);
		end

		iters = iters + 1;
	end

	S = Tr;

end

function [S, cost_sol] = l1averaging(Models, epsilon, maxIters)
	sTr_l2 = l2averaging(Models, epsilon, maxIters);

	%now the L1
	St = sTr_l2;
	done = false;
	iters = 1;
	while(~done && (iters < maxIters) )
		r = zeros(1, 6);
		u = 0;
		sum_norm_vi = 0;
		for i=1:numel(Models)
			Mi = Models{i};
			Tr_i = [Mi; 0 0 0 1];
			vi = se3Log(Tr_i*inv(St));
			norm_vi = norm(vi);

			if(norm_vi > 0)
				sum_norm_vi = sum_norm_vi + norm_vi;
				u = u + 1/norm_vi;
				r = r + vi / norm_vi;
			end
		end

		delta = r / u;
		cost_sol = sum_norm_vi;

		if(norm(delta) < epsilon)
			done = true;
		else
			St = se3Exp(delta)*St;
		end

		iters = iters + 1;

		str = sprintf('L1 cost = %f', sum_norm_vi);
		%disp(str);
	end

	S = St(1:3, 1:4);
end

function indices = selectModels(data, NM)
	[Y, I] = sort(data);

	indices = I(2:NM+1);
end

function [scores] = estimateScores(models, RMM_left, RMM_right)
	N = size(models, 2);
	scores = zeros(N, 1);

	for i=1:N
		scores(i) = evalSolution(models{1, i}, RMM_left, RMM_right);		
	end
end

function residual = evalSolution(Tr, RMM_left, RMM_right)
		X = [Tr(1, 1), Tr(1, 2), Tr(1, 3), Tr(2, 1), Tr(2, 2), Tr(2, 3), Tr(3, 1), Tr(3, 2), Tr(3, 3), Tr(1, 4), Tr(2, 4), Tr(3, 4), 1];
		residual = abs(X*RMM_left*X') +  abs(X*RMM_right*X');
end

function [models, RMM_left, RMM_right, X_left_p, pts_left_c, pts_right_c] = generateModels(Data, K, B, options)
    
	% We are considering the movement from tK+1 to tK
	% but you can change this or just invert the final trnasformation
	pts_left_c = Data(:, 1:2);
	pts_right_c = Data(:, 3:4);
	pts_left_p = Data(:, 5:6);
	pts_right_p = Data(:, 7:8);

	X_left_p = triangulate(pts_left_p, pts_right_p, K, B);
	N = size(X_left_p, 1);

	f = K(1, 1);
	cu = K(1, 3);
	cv = K(2, 3);

	% RMM creation for the trifocal tensor
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	N = size(X_left_p, 1);
	for i=1:N
		X = X_left_p(i, 1);
		Y = X_left_p(i, 2);
		Z = X_left_p(i, 3);

		u1 = pts_left_c(i, 1);
		v1 = pts_left_c(i, 2);

		u2 = pts_right_c(i, 1);
		v2 = pts_right_c(i, 2);

		W_left(3*i-2, :) = [0, 0, 0, f*X, f*Y, f*Z, (cv-v1)*X, (cv-v1)*Y, (cv-v1)*Z, 0, f, cv-v1, 0];
		W_left(3*i-1, :) = [-f*X, -f*Y, -f*Z, 0, 0, 0, (u1-cu)*X, (u1-cu)*Y, (u1-cu)*Z, -f, 0, u1-cu, 0];
		W_left(3*i-0, :) = [f*v1*X, f*v1*Y, f*v1*Z, -f*u1*X, -f*u1*Y, -f*u1*Z, (cu*v1-cv*u1)*X, (cu*v1-cv*u1)*Y, (cu*v1-cv*u1)*Z, f*v1, -f*u1, -cv*u1+cu*v1, 0];

		W_right(3*i-2, :) = [0, 0, 0, f*X, f*Y, f*Z, (cv-v2)*X, (cv-v2)*Y, (cv-v2)*Z, 0, f, cv-v2, 0];
		W_right(3*i-1, :) = [-f*X, -f*Y, -f*Z, 0, 0, 0, (u2-cu)*X, (u2-cu)*Y, (u2-cu)*Z, -f, 0, u2-cu, B*f];
		W_right(3*i-0, :) = [f*v2*X, f*v2*Y, f*v2*Z, -f*u2*X, -f*u2*Y, -f*u2*Z, (cu*v2-cv*u2)*X, (cu*v2-cv*u2)*Y, (cu*v2-cv*u2)*Z, f*v2, -f*u2, -cv*u2+cu*v2, -B*f*v2];
		
	end

	[Ql ,Rl] = qr(W_left, 0);
	RMM_left = Rl;

	[Qr ,Rr] = qr(W_right, 0);
	RMM_right = Rr;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% now we perform a consensus step
	maxiters = options.nModels;
	bestSol = [];
	threshold = 1.5;

	array_sols = cell(maxiters, 1);
	array_res = zeros(maxiters, 1);
	array_in = zeros(maxiters, 1);
    	maxtotalIn = 0;
	iters = 1;
    
	while( iters < maxiters )
		indices = getSamples(N, 3);

		% prepare the optimizer
		myoptions = optimset('Algorithm', {'levenberg-marquardt', 1.0});
		myoptions = optimset(myoptions, 'Jacobian','off');
		myoptions = optimset(myoptions, 'DerivativeCheck', 'off');
		myoptions = optimset(myoptions, 'TolFun', 1e-6);
		myoptions = optimset(myoptions, 'TolX', 1e-9);
		myoptions = optimset(myoptions, 'Display', 'off');
		myoptions = optimset(myoptions, 'MaxFunEvals', 100000);
		myoptions = optimset(myoptions, 'MaxIter', 100);

		Seed = [0, 0, 0, 0, 0, 0.8];
		[solv, resnorm, res, exitflag, output] = lsqnonlin(@costBackproj, Seed, [], [], myoptions, X_left_p(indices, :), pts_left_c(indices, :), pts_right_c(indices, :), K, B);

		R = convert2Cardan(solv(1:3));
		T = solv(4:6)';
		curr_sol = [R, T];
		array_sols{iters, 1} = curr_sol;

		iters = iters + 1;
    	end

	array_sols3 = {};
	idx0 = 1;
	for h=1:numel(array_sols)
		mm = array_sols{h};
		if ~isempty(mm)
			array_sols3{idx0} = mm;
			idx0 = idx0 + 1;
		end
	end
    
  	models = array_sols3;
	
end

function indices = getSamples(N, S)
	indices = randi(N, S, 1);
end

function R = convert2Cardan(v)
	R1 = [1 0 0; 0 cos(v(1)) -sin(v(1)); 0 sin(v(1)) cos(v(1))];
	R2 = [cos(v(2)) 0 -sin(v(2)); 0 1 0; sin(v(2)) 0 cos(v(2))];
	R3 = [cos(v(3)) -sin(v(3)) 0; sin(v(3)) cos(v(3)) 0; 0 0 1];
	
	R = R1*R2*R3;
end

function X = triangulate(pts_left, pts_right, K, B)
	f = K(1, 1);
	cu = K(1, 3);
	cv = K(2, 3);

	X = zeros(size(pts_left, 1), 3);

	for i=1:size(pts_left, 1)
		d = max(pts_left(i, 1) - pts_right(i, 1), 1.0);
		X(i, :) = [(pts_left(i, 1) - cu)*B/d, (pts_left(i, 2) - cv)*B/d, f*B/d ];
	end
end

function [Error, Jac] = costBackproj(v, X_left_p, pts_left_c, pts_right_c, K, B)
		R = convert2Cardan(v(1:3));
		T = v(4:6)';

		Error = zeros(4*size(X_left_p, 1), 1);
		Jac = zeros(4*size(X_left_p, 1), 6);

		% transform the points from the previous left view
		% to the current left view and the current right view
		X_left_c = zeros(size(X_left_p));
		X_right_c = zeros(size(X_left_p));
		for i=1:size(X_left_p, 1)
			X_left_c(i, :) = R*X_left_p(i, :)' + T;
			X_right_c(i, :) = X_left_c(i, :) - [B, 0, 0];
		end

		%now let's compute the error
		JCardan = derivativeCardan(v(1:3));

		for i=1:size(X_left_p, 1)
			expected_lc_u = K(1, 1) * X_left_c(i, 1)/X_left_c(i, 3) + K(1, 3);
			expected_lc_v = K(1, 1) * X_left_c(i, 2)/X_left_c(i, 3) + K(2, 3);
			expected_rc_u = K(1, 1) * X_right_c(i, 1)/X_right_c(i, 3) + K(1, 3);
			expected_rc_v = K(1, 1) * X_right_c(i, 2)/X_right_c(i, 3) + K(2, 3);

			weight = 1.0;
			Error((i*4)-3) = weight*(expected_lc_u - pts_left_c(i, 1));
			Error((i*4)-2) = weight*(expected_lc_v - pts_left_c(i, 2));
			Error((i*4)-1) = weight*(expected_rc_u - pts_right_c(i, 1));
			Error((i*4)-0) = weight*(expected_rc_v - pts_right_c(i, 2));

			for k=1:3
				B = eye(3);
				DR = JCardan{k};
				DR = DR(1:3, 1:3);
				A = DR*[X_left_p(i, :)'];
				Jac((i*4)-3, k) = weight*K(1, 1)*(A(1)*X_left_c(i, 3) - X_left_c(i, 1)*A(3))/(X_left_c(i, 3)*X_left_c(i, 3)); 
     				Jac((i*4)-2, k) = weight*K(1, 1)*(A(2)*X_left_c(i, 3) - X_left_c(i, 2)*A(3))/(X_left_c(i, 3)*X_left_c(i, 3));
     				Jac((i*4)-1, k) = weight*K(1, 1)*(A(1)*X_left_c(i, 3) - X_right_c(i, 1)*A(3))/(X_left_c(i, 3)*X_left_c(i, 3));
     				Jac((i*4)-0, k) = weight*K(1, 1)*(A(2)*X_left_c(i, 3) - X_left_c(i, 2)*A(3))/(X_left_c(i, 3)*X_left_c(i, 3));
			
				Jac((i*4)-3, k+3) = weight*K(1, 1)*(B(k, 1)*X_left_c(i, 3) - X_left_c(i, 1)*B(k, 3))/(X_left_c(i, 3)*X_left_c(i, 3)); 
     				Jac((i*4)-2, k+3) = weight*K(1, 1)*(B(k, 2)*X_left_c(i, 3) - X_left_c(i, 2)*B(k, 3))/(X_left_c(i, 3)*X_left_c(i, 3));
     				Jac((i*4)-1, k+3) = weight*K(1, 1)*(B(k, 1)*X_left_c(i, 3) - X_right_c(i, 1)*B(k, 3))/(X_left_c(i, 3)*X_left_c(i, 3));
     				Jac((i*4)-0, k+3) = weight*K(1, 1)*(B(k, 2)*X_left_c(i, 3) - X_left_c(i, 2)*B(k, 3))/(X_left_c(i, 3)*X_left_c(i, 3));
			end

			
		end
end

function sol = recalculateModelInliers(Tr, X_left_p, pts_left_c, pts_right_c, K, B, threshold)

		R = Tr(1:3, 1:3);
		T = Tr(1:3, 4);
		indicesInliers = false(size(X_left_p, 1), 1);
	
		X_left_c = zeros(size(X_left_p));
		X_right_c = zeros(size(X_left_p));
		for i=1:size(X_left_p, 1)
			X_left_c(i, :) = R*X_left_p(i, :)' + T;
			X_right_c(i, :) = X_left_c(i, :) - [B, 0, 0];
		end

		%now let's compute the error
		inliers = 0;
        
		for i=1:size(X_left_p, 1)
			expected_lc_u = K(1, 1) * X_left_c(i, 1)/X_left_c(i, 3) + K(1, 3);
			expected_lc_v = K(1, 1) * X_left_c(i, 2)/X_left_c(i, 3) + K(2, 3);
			expected_rc_u = K(1, 1) * X_right_c(i, 1)/X_right_c(i, 3) + K(1, 3);
			expected_rc_v = K(1, 1) * X_right_c(i, 2)/X_right_c(i, 3) + K(2, 3);
		
			w = 1.0;
			d = [(pts_left_c(i, 1) - expected_lc_u), (pts_left_c(i, 2) - expected_lc_v), (pts_right_c(i, 1) - expected_rc_u), (pts_right_c(i, 2) - expected_rc_v)];
           		

			if(norm(d) < threshold)
				inliers = inliers + 1;
				indicesInliers(i) = true;
			end
		end

		myoptions = optimset('Algorithm', {'levenberg-marquardt', 1.0});
		myoptions = optimset(myoptions, 'Jacobian','off');
		myoptions = optimset(myoptions, 'DerivativeCheck', 'off');
		myoptions = optimset(myoptions, 'TolFun', 1e-9);
		myoptions = optimset(myoptions, 'TolX', 1e-9);
		myoptions = optimset(myoptions, 'Display', 'off');
		myoptions = optimset(myoptions, 'MaxFunEvals', 100000);
		myoptions = optimset(myoptions, 'MaxIter', 10000);

		Seed = [0, 0, 0, 0, 0, 0];
		[solv, resnorm, res, exitflag, output] = lsqnonlin(@costBackproj, Seed, [], [], myoptions, X_left_p(indicesInliers, :), pts_left_c(indicesInliers, :), pts_right_c(indicesInliers, :), K, B);

		R = convert2Cardan(solv(1:3));
		T = solv(4:6)';
		sol = [R, T];
end
