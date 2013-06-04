#include <mpi.h>
#include "mds.h"
#include "Timer.h"

// 1 - duomenu failas
// 2 - itaraciju skaicius
// 3 - pradiniu vektoriu parinkimo strategija (1 - ant tieses, 2 - atsitiktinis, 3 - didziausiu dispersiju)
// 4 - rezultatu direktorija + failo pradzia (pvz: "/home/Regimantas/tmp/viz_" tai bus sukurta /home/Regmantas/tmp/viz_0.txt, /home/Regmantas/tmp/viz_1.txt ir t.t.)

int main(int argc, char *argv[])
{
    int iMyRank, iProcesors;
    double error;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iMyRank);
    MPI_Comm_size(MPI_COMM_WORLD, &iProcesors);

    //srand((unsigned)time(NULL) * (iMyRank + 1) * 10);

    Mds* pMds = new Mds();
    pMds->SetFile(argv[1]);
    pMds->SetIterations(argv[2]);
    pMds->SetInitMethod(argv[3]);
    pMds->SetRank(iMyRank);
    pMds->SetSeed(((int)time(NULL)) - (iMyRank + 1) * 100);

    time_t tv;
    //Tools::TimerStart(&tv);
	Timer* timer;
	timer = new Timer();
	timer->start();

    pMds->Train(false); //pagal nutylejima ne Seidel

    //double time = Tools::TimerStop(tv);
	double time = timer->stop();

    char *pResultsFile = new char[100];
    sprintf(pResultsFile, "%s", argv[4]);
    char *pRank = new char[4];
    sprintf(pRank, "%d", iMyRank);
    strcat(pResultsFile, pRank);
    strcat(pResultsFile, ".txt");
    pMds->SaveResults(pResultsFile, time);

    delete pMds;
    pMds = NULL;

    MPI_Finalize();

    return 1;
}
