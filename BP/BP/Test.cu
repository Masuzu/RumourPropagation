#include <cuda.h>
#include <cuda_runtime.h>

__global__ void MatrixMulKernel(float * Md, float * Nd, float * Pd, int Width)
{
    // identifiant de thread � deux dimensions, comme la matrice
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    // Pvaleur sert au stockage de la valeur calcul�e par le thread
    float Pvaleur = 0;
    for (int i = 0; i < Width; ++i)
    {
        float MdElement = Md[ty * Width + i];
        float NdElement = Nd[i  * Width + tx];
        Pvaleur        += MdElement * NdElement;
    }
    // �crit la valeur calcul�e dans la matrice de r�sultat
    // chaque thread ne peut �crire qu'une valeur !
    Pd[ty * Width + tx] = Pvaleur;
}

void MatrixMulOnDevice(float * M, float * N, float * P, int Width)
{
    //calcul de la taille des matrices
    int size = Width * Width * sizeof(float);

	float *Md;
	float *Nd;
	float *Pd;

    //allocation des matrices et leur remplissage
    cudaMalloc((void**) &Md, size);
    cudaMemcpy(Md, M, size, cudaMemcpyHostToDevice) ;
    cudaMalloc((void**) &Nd, size);
    cudaMemcpy(Nd, N, size, cudaMemcpyHostToDevice);

    //allocation de la matrice de r�sultat
    cudaMalloc((void**) &Pd, size);

    //multiplication d'une seule matrice
    dim3 dimGrid(1, 1);
    //matrice carr�e
    dim3 dimBlock(Width, Width);

    //produit matriciel proprement dit
    MatrixMulKernel<<<dimGrid, dimBlock>>>(Md, Nd, Pd, Width);

    //r�cup�ration du r�sultat du calcul
    cudaMemcpy(P, Pd, size, cudaMemcpyDeviceToHost);

    //destruction des matrices, d�sormais inutilis�es
    cudaFree(Md);
    cudaFree(Nd);
    cudaFree(Pd);
}