
#include<R.h>
#include<Rmath.h>
#include<R_ext/BLAS.h>
#include<memory.h>

// [rho, resid, sigma] = nchub(A,weight,lambda1,niter_in,niter_out,tol);
char upper = 'U';
char NoTrans = 'N';
char Trans = 'T';
double One = 1.0, MinusOne = -1.0, Zero=0.0;
int IntOne = 1;


double max_val(double a, double b)
{
	double res = a;
	if(a<=b) res = b;
	return res;
}


double sign(double x)
{
	double sgn=0;
	if(x>0) sgn = 1;
	if(x<0) sgn = -1;

	return sgn;
}


// Arguments : X, hub_indx, alpha, lambda1, niter_in, niter_out, tol, output: rho, resid, sigma
void espace(int *N, int *DIM, int *N_HUB,double *X, int *hub_indx, double *ALPHA, double *LAMBDA, int *NITER_IN,int *NITER_OUT, double *TOL,
	           double *rho, double *resid, double *sigma)
{
	double alpha = *ALPHA;
	double *sigma_old;
	double lambda1=*LAMBDA;
	double xe_pij, *x_pii, xe_pji;
	double *temp_coef, temp=0;
	double YX, XX, Xe;
	double max_diff = 0, sig_diff=0;
	int n=*N, p=*DIM;
	int i, j, k, l, q, iter;
	int *indx1, *indx2, length=0;
	int nhub = *N_HUB;
	int *w_indx;

	memset(rho, 0, p*p*sizeof(double));
	memset(resid, 0, n*p*sizeof(double));
	memset(sigma, 0, p*sizeof(double));

	indx1 = (int*)malloc(p*(p-1)/2*sizeof(int));
	indx2 = (int*)malloc(p*(p-1)/2*sizeof(int));
	w_indx = (int*)malloc(p*sizeof(int));
	memset(w_indx, 0, p*sizeof(int));
	for(i=0;i<nhub;i++)
		w_indx[hub_indx[i]-1] = 1;

	sigma_old = (double *)malloc(p*sizeof(double));
    memset(sigma_old, 0, p*sizeof(double));

	memset(resid, 0, n*p*sizeof(double));

	temp_coef = (double *)malloc(p*sizeof(double));
    memset(temp_coef, 0, p*sizeof(double));

	x_pii = (double*)malloc(p*sizeof(double));

	for(i=0;i<p;i++)
	{
		x_pii[i] = F77_NAME(ddot)(&n, &X[n*i], &IntOne, &X[n*i], &IntOne);
		sigma[i] = 1.0;
	}
	memcpy(sigma_old, sigma, p*sizeof(double));

	// Start iteration

	for(iter=0;iter< *NITER_OUT;iter++)
	{

		// Initializatioin of rho

		for(i=0;i<(p-1);i++)
			for(j=i+1;j<p;j++)
			{
				xe_pij = F77_NAME(ddot)(&n, &X[n*i], &IntOne, &X[n*j], &IntOne);

				YX = xe_pij * sqrt(sigma[j]/sigma[i]) + xe_pij*sqrt(sigma[i]/sigma[j]);

				XX = x_pii[j]*sigma[j]/sigma[i] + x_pii[i]*sigma[i]/sigma[j];

				rho[j+p*i] = sign(YX)*(double)max_val((fabs(YX)-lambda1)/XX, 0.0);//(YX/(XX+2*lambda1));

				if(rho[j+p*i]>1) rho[j+p*i]=1;
				if(rho[j+p*i]<-1) rho[j+p*i]=-1;

				rho[i+p*j] = rho[j+p*i];
			}

		// Update residuals

		memcpy(resid, X, n*p*sizeof(double));

		for(i=0;i<p;i++)
		{
			//memcpy(temp_coef,rho_old[i*p],p*sizeof(double));
			for(j=0;j<p;j++)
				temp_coef[j] = rho[j+p*i] * sqrt(sigma[j]/sigma[i]);
			temp_coef[i] = 0;
			F77_NAME(dgemv)(&NoTrans, &n, &p, &MinusOne, X, &n, temp_coef, &IntOne, &One, &resid[n*i], &IntOne);
		}


		for(k=0;k<*NITER_IN;k++)
		{

			//find non-zero differences
			l = 0;
			for(i=0;i<(p-1);i++)
				for(j=i+1;j<p;j++)
				{
					if(rho[i*p+j]!=0)
					{
						indx1[l] = i;
						indx2[l] = j;
						l += 1;
					}
				}

			length = l;

			for(l=0;l<*NITER_IN;l++)
			{
				max_diff = 0;

				for(q = 0;q<length;q++)
				{
					i = (int)indx1[q];
					j = (int)indx2[q];

					//rho_old
					rho[j+p*i] = rho[i+p*j];

					// Y_i^T e_j
					xe_pij = F77_NAME(ddot)(&n, &X[n*i], &IntOne, &resid[n*j], &IntOne);

					// Y_j^T e_i
					xe_pji = F77_NAME(ddot)(&n, &X[n*j], &IntOne, &resid[n*i], &IntOne);

					Xe = xe_pji * sqrt(sigma[j]/sigma[i]) + xe_pij*sqrt(sigma[i]/sigma[j]);

					XX = x_pii[j]*sigma[j]/sigma[i] + x_pii[i]*sigma[i]/sigma[j];

					// beta0
					temp = Xe+rho[j+p*i]*XX;
					if(w_indx[i]==1 || w_indx[j]==1)
						rho[i+p*j] =sign(temp)*(double)max_val((fabs(temp)-lambda1*alpha)/XX, 0.0);// (Xe+rho_old[i+p*j]*XX)/(XX+2*lambda1);
					else
						rho[i+p*j] =sign(temp)*(double)max_val((fabs(temp)-lambda1)/XX, 0.0);// (Xe+rho_old[i+p*j]*XX)/(XX+2*lambda1);

					if(rho[i+p*j]>1) rho[i+p*j]=1;
					if(rho[i+p*j]<-1) rho[i+p*j]=-1;

					
					temp = rho[j+p*i]-rho[i+p*j];
					if(max_diff<fabs(temp))
						max_diff = fabs(temp);

					//Update residual

					// X_i <- X_j
					temp = (rho[j+p*i]-rho[i+p*j])*sqrt(sigma[j]/sigma[i]);
					F77_NAME(daxpy)(&n, &temp, &X[n*j], &IntOne, &resid[n*i], &IntOne);

					// X_j <- X_i
					temp = (rho[j+p*i]-rho[i+p*j])*sqrt(sigma[i]/sigma[j]);
					F77_NAME(daxpy)(&n, &temp, &X[n*i], &IntOne, &resid[n*j], &IntOne);

				}

				if(max_diff<= *TOL)
					break;
			}

			// Shooting on rho
			
			max_diff = 0;
			// update rho
			for(i=0;i<(p-1);i++)
			{
				for(j=i+1;j<p;j++)
				{

					rho[j+p*i] = rho[i+p*j];

					// Y_i^T e_j
					xe_pij = F77_NAME(ddot)(&n, &X[n*i], &IntOne, &resid[n*j], &IntOne);

					// Y_j^T e_i
					xe_pji = F77_NAME(ddot)(&n, &X[n*j], &IntOne, &resid[n*i], &IntOne);

					Xe = xe_pji * sqrt(sigma[j]/sigma[i]) + xe_pij*sqrt(sigma[i]/sigma[j]);

					XX = x_pii[j]*sigma[j]/sigma[i] + x_pii[i]*sigma[i]/sigma[j];

					// beta0
					temp = Xe+rho[j+p*i]*XX;
					if(w_indx[i]==1 || w_indx[j]==1)
						rho[i+p*j] =sign(temp)*(double)max_val((fabs(temp)-lambda1*alpha)/XX, 0.0);// (Xe+rho_old[i+p*j]*XX)/(XX+2*lambda1);
					else
						rho[i+p*j] =sign(temp)*(double)max_val((fabs(temp)-lambda1)/XX, 0.0);// (Xe+rho_old[i+p*j]*XX)/(XX+2*lambda1);

					if(rho[i+p*j]>1) rho[i+p*j]=1;
					if(rho[i+p*j]<-1) rho[i+p*j]=-1;

					
					temp = rho[j+p*i]-rho[i+p*j];
					if(max_diff<fabs(temp))
						max_diff = fabs(temp);

					//Update residual

					// X_i <- X_j
					temp = (rho[j+p*i]-rho[i+p*j])*sqrt(sigma[j]/sigma[i]);
					F77_NAME(daxpy)(&n, &temp, &X[n*j], &IntOne, &resid[n*i], &IntOne);

					// X_j <- X_i
					temp = (rho[j+p*i]-rho[i+p*j])*sqrt(sigma[i]/sigma[j]);
					F77_NAME(daxpy)(&n, &temp, &X[n*i], &IntOne, &resid[n*j], &IntOne);


				}
			}

			if(max_diff<= *TOL)
				break;


		}

		for(i=0;i<(p-1);i++)
			for(j=i+1;j<p;j++)
				rho[j+p*i] = rho[i+p*j];
             
		memcpy(sigma_old, sigma, p*sizeof(double));

		sig_diff = 0;
		for(i=0;i<p;i++)
		{
			sigma[i] =  ((double)n)/ pow(F77_NAME(dnrm2)(&n, &resid[n*i], &IntOne), 2.0);
			temp = fabs(sigma[i]-sigma_old[i]);
			if(sig_diff<temp)
				sig_diff = temp;
		}

		if(sig_diff<= *TOL)
			break;


	}

	//Rprintf("\n Sigma difference = %.4f \n", sig_diff);
	Rprintf("\n All set - Max difference = %.4f ", max_diff);
	Rprintf("\n Iter = %d \n", iter);

	free(w_indx);
	free(temp_coef);
	free(x_pii);
	free(sigma_old);
	free(indx1);
	free(indx2);
}
