#include <mpi.h>
#include "tools.h"
#include "rmds.h"

// 1 - duomenu failas
// 2 - itaraciju skaicius
// 3 - pradiniu vektoriu parinkimo strategija (1 - ant tieses, 2 - atsitiktinis, 3 - didziausiu dispersiju)
// 4 - rezultatu direktorija + failo pradzia (pvz: "/home/Regimantas/tmp/viz_" tai bus sukurta /home/Regmantas/tmp/viz_0.txt, /home/Regmantas/tmp/viz_1.txt ir t.t.)
// 5 - baziniu vektoriu skaicius
int main(int argc, char *argv[])
{
    int iMyRank, iProcesors;
    double error;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iMyRank);
    MPI_Comm_size(MPI_COMM_WORLD, &iProcesors);

    //srand((unsigned)time(NULL) * (iMyRank + 1) * 10);

    Rmds* pRmds = new Rmds();
    pRmds->SetFile(argv[1]); //pRmds->SetFile((char*) "../data3.txt");
    pRmds->SetIterations(argv[2]); //pRmds->SetIterations((char*) "5000");
    pRmds->SetInitMethod(argv[3]); // pRmds->SetInitMethod((char*) "1");
    pRmds->SetRank(iMyRank);
    pRmds->SetSeed(((int)time(NULL)) - (iMyRank + 1) * 100);

    pRmds->SetBasisRows((char*) argv[5]);
    pRmds->SetSizeOfNew((char*) "1");
    pRmds->SetBasisIndextMethod((char*) argv[6]);

    time_t tv;
    Tools::TimerStart(&tv);

    pRmds->Train();

    double time = Tools::TimerStop(tv);

    char *pResultsFile = new char[100];
    sprintf(pResultsFile, "%s", argv[4]);
    char *pRank = new char[4];
    sprintf(pRank, "%d", iMyRank);
    strcat(pResultsFile, pRank);
    strcat(pResultsFile, ".txt");
    pRmds->SaveResults(pResultsFile, time); //pRmds->SaveResults((char*)"/home/Regimantas/tmp/xxx.txt", time, error, pProj);

    delete pRmds; //free
    pRmds = NULL;

    MPI_Finalize();

    return 1;
}
