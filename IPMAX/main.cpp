//
///
/// Created by DBL-XQ on 2024/10/15
//

#include "stdafx.h"
#include "argument.h"
#include <omp.h>

#include "graph.h"
#include "TemporalGraph.h"
#include "InfPerMAX.h"
#include "compare.h"
#include "sccdag_gen.h"

void OutputSeedSetToFile(vector<int> seed_set)
{
    ofstream of;
    of.open("seedfile");
    for (int seed: seed_set)
    {
        of << seed << " ";
    }
    of << endl;
    of.close();
}


void run_with_parameter(TemporalGraph & tg, const Argument & arg){



    if (arg._alg == "RG"){
        double startTime = omp_get_wtime();
        IPMAX::ReverseGreedy(tg, arg);
        double endTime = omp_get_wtime();
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;
        cout << " R: " << arg._R << " W: " << arg._W << " k: " << arg._k << " theta: " << arg._theta << " alpha: " << arg._alpha << " threshold nodes: " << arg._alpha * tg.V<< endl;
        cout << "******************************** ReverseGreedy end ***************************************" << endl;
    }


    else if (arg._alg == "LR"){
        double startTime = omp_get_wtime();
        IPMAX::LazyReplace_prn_opt(tg, arg);
        double endTime = omp_get_wtime();
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;
        cout << " R: " << arg._R << " W: " << arg._W << " k: " << arg._k << " theta: " << arg._theta << " alpha: " << arg._alpha << " threshold nodes: " << arg._alpha * tg.V<< endl;
        cout << "******************************** LazyReplace end ***************************************" << endl;
    }

    else if (arg._alg == "RGSCC"){
        unordered_map<int, unordered_map<int,vector<int>>> FT_sccs;                         // [t][scc_id][node]
        unordered_map<int, unordered_map<int, vector<pair<int, double>>>> fin_all_scc_dags;  // [t][scc_id] -> list of edges
        SCC sccdag;
        sccdag.generate_scc_dags(tg, arg, FT_sccs, fin_all_scc_dags);
        double startTime = omp_get_wtime();
        IPMAX::ReverseGreedySCC(tg, arg,  FT_sccs, fin_all_scc_dags);
        double endTime = omp_get_wtime();
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;
        cout << " R: " << arg._R << " W: " << arg._W << " k: " << arg._k << " theta: " << arg._theta << " alpha: " << arg._alpha << " threshold nodes: " << arg._alpha * tg.V<< endl;
        cout << "******************************** ReverseGreedy-SCC end ***************************************" << endl;
    }


    else if (arg._alg == "LRSCC"){
        unordered_map<int, unordered_map<int,vector<int>>> FT_sccs;                        // [t][scc_id][node]
        unordered_map<int, unordered_map<int, vector<pair<int, double>>>> fin_all_scc_dags;  // [t][scc_id] -> list of edges
        SCC sccdag;
        sccdag.generate_scc_dags(tg, arg, FT_sccs,  fin_all_scc_dags);
        double startTime = omp_get_wtime();
        IPMAX::LazyReplaceSCC(tg, arg,  FT_sccs, fin_all_scc_dags);
        double endTime = omp_get_wtime();
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;
        cout << " R: " << arg._R << " W: " << arg._W << " k: " << arg._k << " theta: " << arg._theta << " alpha: " << arg._alpha << " threshold nodes: " << arg._alpha * tg.V<< endl;
        cout << "******************************** LazyReplace-SCC end ***************************************" << endl;
    }

    // in window W continuous

    else if (arg._alg == "WRG"){
        double startTime = omp_get_wtime();
        IPMAX::WinReverseGreedy(tg, arg);
        double endTime = omp_get_wtime();
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;
        cout << " R: " << arg._R << " W: " << arg._W << " k: " << arg._k << " theta: " << arg._theta << " alpha: " << arg._alpha << " threshold nodes: " << arg._alpha * tg.V<< endl;
        cout << "******************************** WinReverseGreedy end ***************************************" << endl;
    }

    else if (arg._alg == "WLR"){
        double startTime = omp_get_wtime();
        IPMAX::WinLazyReplace(tg, arg);
        double endTime = omp_get_wtime();
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;
        cout << " R: " << arg._R << " W: " << arg._W << " k: " << arg._k << " theta: " << arg._theta << " alpha: " << arg._alpha << " threshold nodes: " << arg._alpha * tg.V<< endl;
        cout << "******************************** WinLazyReplace end ***************************************" << endl;
    }

    else{
        cout << "wrong algorithm name." <<endl;
    }

}


int main(int argc, char * argv[]) {


    if (argc == 1){
        cout << "Usage: ./IPMAX -dataset=../dataset/test/ -k=2 -theta=2 -alpha=0.5" << endl;
        exit(0);
    }

    const Argument arg(argc, argv);


    string attribute_file, graph_file, target_file;
    attribute_file = arg._dataset + "attribute.txt";
    graph_file = arg._dataset + "graph.txt";



    TemporalGraph tg(attribute_file, graph_file);

//    if (int(tg.V * arg._k) > arg._alpha * tg.V){
//        cout << "k (the number of seed nodes is too large)" << endl;
//        exit(0);
//    }

    run_with_parameter(tg, arg);

    return 0;
}
