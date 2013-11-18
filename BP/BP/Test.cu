#include <cuda.h>
#include <cuda_runtime.h>

__global__ void MatrixMulKernel(float * Md, float * Nd, float * Pd, int Width)
{
    // identifiant de thread à deux dimensions, comme la matrice
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    // Pvaleur sert au stockage de la valeur calculée par le thread
    float Pvaleur = 0;
    for (int i = 0; i < Width; ++i)
    {
        float MdElement = Md[ty * Width + i];
        float NdElement = Nd[i  * Width + tx];
        Pvaleur        += MdElement * NdElement;
    }
    // écrit la valeur calculée dans la matrice de résultat
    // chaque thread ne peut écrire qu'une valeur !
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

    //allocation de la matrice de résultat
    cudaMalloc((void**) &Pd, size);

    //multiplication d'une seule matrice
    dim3 dimGrid(1, 1);
    //matrice carrée
    dim3 dimBlock(Width, Width);

    //produit matriciel proprement dit
    MatrixMulKernel<<<dimGrid, dimBlock>>>(Md, Nd, Pd, Width);

    //récupération du résultat du calcul
    cudaMemcpy(P, Pd, size, cudaMemcpyDeviceToHost);

    //destruction des matrices, désormais inutilisées
    cudaFree(Md);
    cudaFree(Nd);
    cudaFree(Pd);
}