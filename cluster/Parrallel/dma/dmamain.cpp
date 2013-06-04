#include <mpi.h>

#include "dma.h"

// 1 - duomenu failas
// 2 - itaraciju skaicius
// 3 - pradiniu vektoriu parinkimo strategija (1 - ant tieses, 2 - atsitiktinis, 3 - didziausiu dispersiju)
// 4 - rezultatu direktorija + failo pradzia (pvz: "/home/Regimantas/tmp/viz_" tai bus sukurta /home/Regmantas/tmp/viz_0.txt, /home/Regmantas/tmp/viz_1.txt ir t.t.)
// 5 - parametras k
int main(int argc, char *argv[])
{
    int iMyRank = 4, iProcesors;
    double error;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iMyRank);
    MPI_Comm_size(MPI_COMM_WORLD, &iProcesors);

    //srand((unsigned)time(NULL) * (iMyRank + 1) * 10);

    Dma *pDma = new Dma();
    pDma->SetFile(argv[1]);
    pDma->SetIterations(argv[2]);
    pDma->SetInitMethod(argv[3]);
    pDma->SetNeighbour(argv[5]);
    pDma->SetRank(iMyRank);
    pDma->SetSeed(((int) time(NULL)) - (iMyRank + 1) * 10);

    time_t tv;
    Tools::TimerStart(&tv);

    pDma->Train(true);

    double time = Tools::TimerStop(tv);

    char *pResultsFile = new char[100];
   
    sprintf(pResultsFile, "%s", argv[4]);
    char *pRank = new char[4];
    sprintf(pRank, "%d", iMyRank);
    strcat(pResultsFile, pRank);
    strcat(pResultsFile, ".txt");

    pDma->SaveResults(pResultsFile, time);

    delete pDma;
    pDma = NULL;

    MPI_Finalize();

    return 1;
}
