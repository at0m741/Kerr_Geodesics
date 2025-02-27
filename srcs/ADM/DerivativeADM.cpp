#include <Geodesics.h>

void fill_ghost_zones() {
    for (int g = 0; g < GHOST; g++) {
        for (int j = 0; j < NY_TOTAL; j++) {
            for (int k = 0; k < NZ_TOTAL; k++) {
                globalGrid[g][j][k] = globalGrid[NX_TOTAL - 2*GHOST + g][j][k];
                globalGrid[NX_TOTAL - GHOST + g][j][k] = globalGrid[GHOST + g][j][k];
            }
        }
    }
}






 double partialX_alpha(int i, int j, int k) {
	if(i < 2 || i > NX-3) {
		return (globalGrid[i+1][j][k].alpha - globalGrid[i-1][j][k].alpha) / (2.0 * DX);
	}
	double alpha_ip2 = globalGrid[i+2][j][k].alpha;
	double alpha_ip1 = globalGrid[i+1][j][k].alpha;
	double alpha_im1 = globalGrid[i-1][j][k].alpha;
	double alpha_im2 = globalGrid[i-2][j][k].alpha;

	return ( -alpha_ip2 + 8.0*alpha_ip1 - 8.0*alpha_im1 + alpha_im2 ) / (12.0 * DX );
}


 double partialY_alpha(int i, int j, int k) {
	if (j < 2 || j > NY-3) {
		return (globalGrid[i][j+1][k].alpha - globalGrid[i][j-1][k].alpha) / (2.0 * DY);
	} 
	double alpha_jp2 = globalGrid[i][j+2][k].alpha;
	double alpha_jp1 = globalGrid[i][j+1][k].alpha;
	double alpha_jm1 = globalGrid[i][j-1][k].alpha;
	double alpha_jm2 = globalGrid[i][j-2][k].alpha;

	return ( -alpha_jp2 + 8.0*alpha_jp1 - 8.0*alpha_jm1 + alpha_jm2 ) / (12.0 * DY );
}


 double partialZ_alpha(int i, int j, int k) {
	if (k < 2 || k > NZ-3) {
		return (globalGrid[i][j][k+1].alpha - globalGrid[i][j][k-1].alpha) / (2.0 * DZ);
	}
	double alpha_kp2 = globalGrid[i][j][k+2].alpha;
	double alpha_kp1 = globalGrid[i][j][k+1].alpha;
	double alpha_km1 = globalGrid[i][j][k-1].alpha;
	double alpha_km2 = globalGrid[i][j][k-2].alpha;

	return ( -alpha_kp2 + 8.0*alpha_kp1 - 8.0*alpha_km1 + alpha_km2 ) / (12.0 * DZ );
}


 double partialXX_alpha(int i, int j, int k) {
	return (globalGrid[i+1][j][k].alpha - 2.0*globalGrid[i][j][k].alpha + globalGrid[i-1][j][k].alpha)/(DX-1*DX-1);
}


 double partialYY_alpha(int i, int j, int k) {
	return (globalGrid[i][j+1][k].alpha - 2.0*globalGrid[i][j][k].alpha + globalGrid[i][j-1][k].alpha)/(DY-1*DY-1);
}


 double partialZZ_alpha(int i, int j, int k) {
	return (globalGrid[i][j][k+1].alpha - 2.0*globalGrid[i][j][k].alpha + globalGrid[i][j][k-1].alpha)/(DZ-1*DZ-1);
}


double partialX_gamma(int i, int j, int k, int a, int b)
{
    if (i >= 2 && i <= NX - 3)
    {
        double res = fourth_order_diff(
            globalGrid[i+2][j][k].gamma[a][b],
            globalGrid[i+1][j][k].gamma[a][b],
            globalGrid[i-1][j][k].gamma[a][b],
            globalGrid[i-2][j][k].gamma[a][b],
            DX 
        );
		return res;
    }
    else if (i >= 1 && i <= NX - 2)
    {
        double res = second_order_diff(
            globalGrid[i+1][j][k].gamma[a][b],
            globalGrid[i-1][j][k].gamma[a][b],
            DX 
        );
        return res;
    }
    else if (i == 0) 
    {
        return (globalGrid[i+1][j][k].gamma[a][b] - globalGrid[i][j][k].gamma[a][b]) / DX;
    }
    else if (i == NX - 1) 
    {
        return (globalGrid[i][j][k].gamma[a][b] - globalGrid[i-1][j][k].gamma[a][b]) / DX;
    }
    return 0.0;
}


double partialY_gamma(int i, int j, int k, int a, int b)
{
	if (j >= 2 && j <= NY - 3)
	{
		return fourth_order_diff(
			globalGrid[i][j+2][k].gamma[a][b],
			globalGrid[i][j+1][k].gamma[a][b],
			globalGrid[i][j-1][k].gamma[a][b],
			globalGrid[i][j-2][k].gamma[a][b],
			DY
		);
	}
	else if (j >= 1 && j <= NY - 2)
	{
		return second_order_diff(
			globalGrid[i][j+1][k].gamma[a][b],
			globalGrid[i][j-1][k].gamma[a][b],
			DY
		);
	}
	else if (j == 0)
	{
		return (globalGrid[i][j+1][k].gamma[a][b] - globalGrid[i][j][k].gamma[a][b]) / DY;
	}
	else if (j == NY - 1)
	{
		return (globalGrid[i][j][k].gamma[a][b] - globalGrid[i][j-1][k].gamma[a][b]) / DY;
	}
	return 0.0;
}

double partialZ_gamma(int i, int j, int k, int a, int b)
{
    if (k >= 2 && k <= NZ - 3)
    {
        return fourth_order_diff(
            globalGrid[i][j][k+2].gamma[a][b],
            globalGrid[i][j][k+1].gamma[a][b],
            globalGrid[i][j][k-1].gamma[a][b],
            globalGrid[i][j][k-2].gamma[a][b],
            DZ
        );
    }
    else if (k >= 1 && k <= NZ - 2)
    {
        return second_order_diff(
            globalGrid[i][j][k+1].gamma[a][b],
            globalGrid[i][j][k-1].gamma[a][b],
            DZ
        );
    }
    else if (k == 0)
	{
		return (globalGrid[i][j][k+1].gamma[a][b] - globalGrid[i][j][k].gamma[a][b]) / DZ;
	}
	else if (k == NZ - 1)
	{
		return (globalGrid[i][j][k].gamma[a][b] - globalGrid[i][j][k-1].gamma[a][b]) / DZ;
	}
	return 0.0;
}



 double partialX_betacomp(int i, int j, int k, int comp) {
	 return (globalGrid[i+1][j][k].beta[comp] 
			 - globalGrid[i-1][j][k].beta[comp]) / (2.0 * DX);
}


 double partialY_betacomp(int i, int j, int k, int comp) {
	return (globalGrid[i][j+1][k].beta[comp] 
			- globalGrid[i][j-1][k].beta[comp]) / (2.0 * DY);
}


 double partialZ_betacomp(int i, int j, int k, int comp) {
	return (globalGrid[i][j][k+1].beta[comp] 
			- globalGrid[i][j][k-1].beta[comp])/(2.0 * DZ);
}


double partialXY_alpha(int i, int j, int k) {
    return ( globalGrid[i+1][j+1][k].alpha - globalGrid[i+1][j-1][k].alpha
           - globalGrid[i-1][j+1][k].alpha + globalGrid[i-1][j-1][k].alpha )
           / (4.0 * DX * DY);
}

double partialXZ_alpha(int i, int j, int k) {
	return ( globalGrid[i+1][j][k+1].alpha - globalGrid[i+1][j][k-1].alpha
		   - globalGrid[i-1][j][k+1].alpha + globalGrid[i-1][j][k-1].alpha )
		   / (4.0 * DX * DZ);
}

double partialYZ_alpha(int i, int j, int k) {
	return ( globalGrid[i][j+1][k+1].alpha - globalGrid[i][j+1][k-1].alpha
		   - globalGrid[i][j-1][k+1].alpha + globalGrid[i][j-1][k-1].alpha )
		   / (4.0 * DY * DZ);
}

double second_partial_alpha(int i, int j, int k, int a, int b)
{
    if(a==0 && b==0) return partialXX_alpha(i,j,k);
    if(a==1 && b==1) return partialYY_alpha(i,j,k);
    if(a==2 && b==2) return partialZZ_alpha(i,j,k);

    if((a==0 && b==1)||(a==1 && b==0)) return partialXY_alpha(i,j,k);
    if((a==0 && b==2)||(a==2 && b==0)) return partialXZ_alpha(i,j,k);
    if((a==1 && b==2)||(a==2 && b==1)) return partialYZ_alpha(i,j,k);

    return 0.0;
}

