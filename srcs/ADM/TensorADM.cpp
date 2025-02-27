#include <Geodesics.h>

static void print_matrix_2D(const char* name, const double mat[3][3])
{
    printf("\n%s (3x3):\n", name);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%12.6f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}


double GridTensor::partialX_gamma(int i, int j, int k, int a, int b)
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


double GridTensor::partialY_gamma(int i, int j, int k, int a, int b)
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

double GridTensor::partialZ_gamma(int i, int j, int k, int a, int b)
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


void GridTensor::compute_partial_christoffel(int i, int j, int k, int dim, double partialGamma[3][3][3][3], double d) {
double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];
    int a = (dim == 0) ? i : ((dim == 1) ? j : k);
    int max_a = (dim == 0) ? NX : ((dim == 1) ? NY : NZ);

    if (a >= 2 && a <= max_a - 3) {
        compute_christoffel_3D(i - 2 * (dim == 0), j - 2 * (dim == 1), k - 2 * (dim == 2), Gmm);
        compute_christoffel_3D(i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
        compute_christoffel_3D(i + 2 * (dim == 0), j + 2 * (dim == 1), k + 2 * (dim == 2), Gpp);
    } else if (a >= 1 && a <= max_a - 2) {
        compute_christoffel_3D(i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
    }

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                if (a >= 2 && a <= max_a - 3) {
                    partialGamma[dim][kk][aa][bb] = fourth_order_diff(Gpp[kk][aa][bb], Gp[kk][aa][bb], Gm[kk][aa][bb], Gmm[kk][aa][bb], d);
                } else if (a >= 1 && a <= max_a - 2) {
                    partialGamma[dim][kk][aa][bb] = second_order_diff(Gp[kk][aa][bb], Gm[kk][aa][bb], d);
                } else {
                    partialGamma[dim][kk][aa][bb] = 0.0;
                }
            }
        }
    }
}

void Grid::compute_ricci_3D(int i, int j, int k, double Ricci[3][3]) {
    double Gamma[3][3][3];
	GridTensor gridTensor;
    gridTensor.compute_christoffel_3D(i, j, k, Gamma);
    double partialGamma[3][3][3][3] = {};

    gridTensor.compute_partial_christoffel(i, j, k, 0, partialGamma, DX);
    gridTensor.compute_partial_christoffel(i, j, k, 1, partialGamma, DY);
    gridTensor.compute_partial_christoffel(i, j, k, 2, partialGamma, DZ);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            for (int m = 0; m < 3; m++) {
                term1 += partialGamma[m][m][a][b];
                term2 += partialGamma[a][m][m][b];
            }
            for (int kk = 0; kk < 3; kk++) {
                for (int ll = 0; ll < 3; ll++) {
                    term3 += Gamma[kk][a][b] * Gamma[ll][kk][ll];
                    term4 += Gamma[ll][a][kk] * Gamma[kk][b][ll];
                }
            }
            Ricci[a][b] = term1 - term2 + term3 - term4;
        }
    }
	/* print_matrix_2D("Ricci", Ricci); */
}
