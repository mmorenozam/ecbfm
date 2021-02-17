functions {
	real[] 	cbf(real t,
				real[] x,
				real[] theta,
				real[] y_r,
				int[] y_i) {
		real dxdt[9];
		real v1;
		real v2;
		real v3;
		real v4;
		real v5;
		real v6;
		real v7;
		real v8;
		
		v1 = theta[12]*exp(-theta[25]/x[9])*x[1]*x[6]/(theta[17]+x[1]);
		v2 = theta[13]*exp(-theta[25]/x[9])*x[2]*x[6]/(theta[18]+x[2]);
		v3 = theta[14]*exp(-theta[25]/x[9])*x[1]*x[7]/(theta[19]+x[1]);
		v4 = theta[15]*exp(-theta[25]/x[9])*x[3]*x[8]/(theta[20]+x[3]);
		v5 = theta[16]*exp(-theta[25]/x[9])*x[4]*x[8]/(theta[21]*x[8]+x[4]);
		v6 = theta[22]*exp(-theta[26]/x[9])*x[6]*x[3];
		v7 = theta[23]*exp(-theta[26]/x[9])*x[7]*x[4];
		v8 = theta[24]*exp(-theta[26]/x[9])*x[8]*x[5]^2;
		
		
		dxdt[1] = - theta[1]*v1 - theta[2]*v3;
		dxdt[2] = - theta[3]*v2;
		dxdt[3] = theta[4]*v1 + theta[5]*v2 - theta[6]*v4;
		dxdt[4] = theta[7]*v3 - theta[8]*v5;
		dxdt[5] = theta[9]*v3 + theta[10]*v4 +theta[11]*v5;
		dxdt[6] = v1 + v2 - v6;
		dxdt[7] = v3 - v7;
		dxdt[8] = v4 + v5 - v8;
		dxdt[9] = theta[29]*(theta[1]*v1+theta[2]*v3) + theta[30]*theta[3]*v2 + theta[31]*theta[6]*v4 + theta[32]*theta[8]*v5 - theta[27]*(x[9] - theta[28]);
		
		return dxdt;
	}
}

data {
	int<lower=1> T;
	real<lower=0> x[T,9];
	real t0;
	real ts[T];
	real x0[9];
	real scl[9];
	real tmax;

}

transformed data {
	real y_r[0];
	int y_i[0];
	real<lower=0> x0_1[9];
	real<lower=0> xn[T,9];
	for (t in 1:T)
		for (n in 1:9)
			xn[t,n]=x[t,n]/scl[n];
  
  for (n in 1:9)
    x0_1[n] = x0[n]/scl[n]; 
}

parameters {
	real <lower=0> sigma;
	real <lower=0> mu1;
	real <lower=0> mu2;
	real <lower=0> mu3;
	real <lower=0> mu4;
	real <lower=0> mu5;
	real <lower=0> ks1;
	real <lower=0> ks2;
	real <lower=0> ks3;
	real <lower=0> ks4;
	real <lower=0> ks5;
	real <lower=0> k1;
	real <lower=0> k2;
	real <lower=0> k3;
	real <lower=0> yc1;
	real <lower=0> yc2;
	real <lower=0> yc3;
	real <lower=0> yc4;
	real <lower=0> yc5;
	real <lower=0> yc6;
	real <lower=0> yc7;
	real <lower=0> yc8; 
	real <lower=0> yc9;
	real <lower=0> yc10;
	real <lower=0> yc11;
	real <lower=0> a;
	real <lower=0> b;
	real <lower=0> ql;
	real <lower=0> te;
	real <lower=0> yq1;
	real <lower=0> yq2;
	real <lower=0> yq3;
	real <lower=0> yq4;
}

transformed parameters {
	real x_hat[T,9];
	{
		real theta[32];
		theta[1] = yc1;
		theta[2] = yc2;
		theta[3] = yc3;
		theta[4] = yc4;
		theta[5] = yc5;
		theta[6] = yc6;
		theta[7] = yc7;
		theta[8] = yc8;
		theta[9] = yc9;
		theta[10] = yc10;
		theta[11] = yc11;
		theta[12] = mu1;
		theta[13] = mu2;
		theta[14] = mu3;
		theta[15] = mu4;
		theta[16] = mu5;
		theta[17] = ks1;
		theta[18] = ks2;
		theta[19] = ks3;
		theta[20] = ks4;
		theta[21] = ks5;
		theta[22] = k1;
		theta[23] = k2;
		theta[24] = k3;
		theta[25] = a;
		theta[26] = b;
		theta[27] = ql;
		theta[28] = te;
		theta[29] = yq1;
		theta[30] = yq2;
		theta[31] = yq3;
		theta[32] = yq4;

		x_hat = integrate_ode_rk45(cbf, x0_1, t0, ts, theta, y_r, y_i,1.0E-6, 1.0E-6, 1.0E4);
	}
}
model{
	mu1~normal(0.5,0.3);
	mu2~normal(0.5,0.3);
	mu3~normal(0.5,0.3);
	mu4~normal(0.5,0.3);
	mu5~normal(0.5,0.3);
	ks1~normal(0.5,0.3);
	ks2~normal(0.5,0.3);
	ks3~normal(0.5,0.3);
	ks4~normal(0.5,0.3);
	ks5~normal(0.5,0.3);
	k1~normal(0.5,0.3);
	k2~normal(0.5,0.3);
	k3~normal(0.5,0.3);
	yc1~normal(0.5,0.3);
	yc2~normal(0.5,0.3);
	yc3~normal(0.5,0.3);
	yc4~normal(0.5,0.3);
	yc5~normal(0.5,0.3);
	yc6~normal(0.5,0.3);
	yc7~normal(0.5,0.3);
	yc8~normal(0.5,0.3);
	yc9~normal(0.5,0.3);
	yc10~normal(0.5,0.3);
	yc11~normal(0.5,0.3);
	a~normal(0.5,0.3);
	b~normal(0.5,0.3);
	ql~normal(0.5,0.3);
	te~normal(tmax,0.01);
	yq1~normal(0.5,0.3);
	yq2~normal(0.5,0.3);
	yq3~normal(0.5,0.3);
	yq4~normal(0.5,0.3);
	sigma~cauchy(0,1);
	for (t in 1:T)
		xn[t]~normal(x_hat[t],sigma);
}


generated quantities{
  real log_lik[T,8];
  for (t in 1:T)
    for (n in 1:8)
      log_lik[t,n] = normal_lpdf(xn[t,n]|x_hat[t,n],sigma);
}

