#include <Geodesics.h>


double partialX_alpha(Grid &grid_obj, int i, int j, int k) {
    if (i < 2 || i > NX - 3) {
        return (grid_obj.getCell(i + 1, j, k).alpha - grid_obj.getCell(i - 1, j, k).alpha) / (2.0 * DX);
    }
    double alpha_ip2 = grid_obj.getCell(i + 2, j, k).alpha;
    double alpha_ip1 = grid_obj.getCell(i + 1, j, k).alpha;
    double alpha_im1 = grid_obj.getCell(i - 1, j, k).alpha;
    double alpha_im2 = grid_obj.getCell(i - 2, j, k).alpha;

    return (-alpha_ip2 + 8.0 * alpha_ip1 - 8.0 * alpha_im1 + alpha_im2) / (12.0 * DX);
}

double partialY_alpha(Grid &grid_obj, int i, int j, int k) {
    if (j < 2 || j > NY - 3) {
        return (grid_obj.getCell(i, j + 1, k).alpha - grid_obj.getCell(i, j - 1, k).alpha) / (2.0 * DY);
    }
    double alpha_jp2 = grid_obj.getCell(i, j + 2, k).alpha;
    double alpha_jp1 = grid_obj.getCell(i, j + 1, k).alpha;
    double alpha_jm1 = grid_obj.getCell(i, j - 1, k).alpha;
    double alpha_jm2 = grid_obj.getCell(i, j - 2, k).alpha;

    return (-alpha_jp2 + 8.0 * alpha_jp1 - 8.0 * alpha_jm1 + alpha_jm2) / (12.0 * DY);
}

double partialZ_alpha(Grid &grid_obj, int i, int j, int k) {
    if (k < 2 || k > NZ - 3) {
        return (grid_obj.getCell(i, j, k + 1).alpha - grid_obj.getCell(i, j, k - 1).alpha) / (2.0 * DZ);
    }
    double alpha_kp2 = grid_obj.getCell(i, j, k + 2).alpha;
    double alpha_kp1 = grid_obj.getCell(i, j, k + 1).alpha;
    double alpha_km1 = grid_obj.getCell(i, j, k - 1).alpha;
    double alpha_km2 = grid_obj.getCell(i, j, k - 2).alpha;

    return (-alpha_kp2 + 8.0 * alpha_kp1 - 8.0 * alpha_km1 + alpha_km2) / (12.0 * DZ);
}


double partialXX_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i+1, j, k).alpha - 2.0 * grid_obj.getCell(i, j, k).alpha + grid_obj.getCell(i-1, j, k).alpha) 
           / (DX * DX);
}

double partialYY_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i, j+1, k).alpha - 2.0 * grid_obj.getCell(i, j, k).alpha + grid_obj.getCell(i, j-1, k).alpha) 
           / (DY * DY);
}

double partialZZ_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i, j, k+1).alpha - 2.0 * grid_obj.getCell(i, j, k).alpha + grid_obj.getCell(i, j, k-1).alpha) 
           / (DZ * DZ);
}


double partialX_betacomp(Grid &grid_obj, int i, int j, int k, int comp) {
    return (grid_obj.getCell(i+1, j, k).beta[comp] - grid_obj.getCell(i-1, j, k).beta[comp]) / (2.0 * DX);
}

double partialY_betacomp(Grid &grid_obj, int i, int j, int k, int comp) {
    return (grid_obj.getCell(i, j+1, k).beta[comp] - grid_obj.getCell(i, j-1, k).beta[comp]) / (2.0 * DY);
}

double partialZ_betacomp(Grid &grid_obj, int i, int j, int k, int comp) {
    return (grid_obj.getCell(i, j, k+1).beta[comp] - grid_obj.getCell(i, j, k-1).beta[comp]) / (2.0 * DZ);
}

double partialXY_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i+1, j+1, k).alpha - grid_obj.getCell(i+1, j-1, k).alpha
            - grid_obj.getCell(i-1, j+1, k).alpha + grid_obj.getCell(i-1, j-1, k).alpha)
           / (4.0 * DX * DY);
}

double partialXZ_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i+1, j, k+1).alpha - grid_obj.getCell(i+1, j, k-1).alpha
            - grid_obj.getCell(i-1, j, k+1).alpha + grid_obj.getCell(i-1, j, k-1).alpha)
           / (4.0 * DX * DZ);
}

double partialYZ_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i, j+1, k+1).alpha - grid_obj.getCell(i, j+1, k-1).alpha
            - grid_obj.getCell(i, j-1, k+1).alpha + grid_obj.getCell(i, j-1, k-1).alpha)
           / (4.0 * DY * DZ);
}

double second_partial_alpha(Grid &grid_obj, int i, int j, int k, int a, int b)
{
	if(a==0 && b==0) return partialXX_alpha(grid_obj, i, j, k);
	if(a==1 && b==1) return partialYY_alpha(grid_obj, i, j, k);
	if(a==2 && b==2) return partialZZ_alpha(grid_obj, i, j, k);

	if((a==0 && b==1)||(a==1 && b==0)) return partialXY_alpha(grid_obj, i, j, k);
	if((a==0 && b==2)||(a==2 && b==0)) return partialXZ_alpha(grid_obj, i, j, k);
	if((a==1 && b==2)||(a==2 && b==1)) return partialYZ_alpha(grid_obj, i, j, k);

	return 0.0;
}

