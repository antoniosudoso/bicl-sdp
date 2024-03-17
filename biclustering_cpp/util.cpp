#include "util.h"
#include "Node.h"
#include "config_params.h"
#include <iomanip>

void save_X_to_file(arma::sp_mat &Xu, arma::sp_mat &Xv){

    std::ofstream f;
    f.open(result_path);
    for (size_t i = 0; i < Xu.n_rows; i++){
        for (size_t j = 0; j < Xu.n_cols; j++){
            if (Xu(i, j) > 0)
                f << " " << 1;
            else
                f << " " << 0;
        }
        f << "\n";
    }
    f << "\n";
    for (size_t i = 0; i < Xv.n_rows; i++){
        for (size_t j = 0; j < Xv.n_cols; j++){
            if (Xv(i, j) > 0)
                f << " " << 1;
            else
                f << " " << 0;
        }
        f << "\n";
    }
    f.close();
}


void print_header_sdp() {

    log_file << "\n" << "|" <<
              std::setw(8) << "N" << "|" <<
              std::setw(8) << "M" << "|" <<
              std::setw(8) << "ID_PAR" << "|" <<
              std::setw(8) << "ID" << "|" <<
              std::setw(12) << "UB_PAR" << "|" <<
              std::setw(12) << "UB" << "|" <<
              std::setw(10) << "TIME (s)" << "|" <<
              std::setw(8) << "CP_ITER" << "|" <<
              std::setw(8) << "CP_FLAG" << "|" <<
              std::setw(10) << "CP_INEQ" << "|" <<
              std::setw(12) << "LB" << "|" <<
              std::setw(12) << "BEST_LB" << "|" <<
              std::setw(8) << "SET" << "|" <<
              std::setw(8) << "I" << " " <<
              std::setw(8) << "J" << "|" <<
              std::setw(13) << "NODE_GAP" << "|" <<
              std::setw(13) << "GAP" << "|" <<
              std::setw(8) << "OPEN" << "|"
              << std::endl;

}

void print_log_sdp(int n, int m, int node_parent, int node, double ub_parent, double ub, double time,
                   int cp_iter, int cp_flag, int n_ineq, double lb, double glb, int type, int i, int j,
                   double node_gap, double gap, int open, bool update) {

    std::string set;
    switch (type) {
        case BRANCH_ON_U:
            set = "U";
            break;
        case BRANCH_ON_V:
            set = "V";
            break;
        default:
            set = "-1";
            break;
    }

    if (!update) {

        log_file << "|" <<
                  std::setw(8) << n << "|" <<
                  std::setw(8) << m << "|" <<
                  std::setw(8) << node_parent << "|" <<
                  std::setw(8) << node << "|" <<
                  std::setw(12) << ub_parent << "|" <<
                  std::setw(12) << ub << "|" <<
                  std::setw(10) << time << "|" <<
                  std::setw(8) << cp_iter << "|" <<
                  std::setw(8) << cp_flag << "|" <<
                  std::setw(10) << n_ineq << "|" <<
                  std::setw(12) << lb << "|" <<
                  std::setw(12) << glb << "|" <<
                  std::setw(8) << set << "|" <<
                  std::setw(8) << i << " " <<
                  std::setw(8) << j << "|" <<
                  std::setw(13) << node_gap << "|" <<
                  std::setw(13) << gap << "|" <<
                  std::setw(8) << open << "|"
                  << std::endl;

    } else {

        log_file << "|" <<
                  std::setw(8) << n << "|" <<
                  std::setw(8) << m << "|" <<
                  std::setw(8) << node_parent << "|" <<
                  std::setw(8) << node << "|" <<
                  std::setw(12) << ub_parent << "|" <<
                  std::setw(12) << ub << "|" <<
                  std::setw(10) << time << "|" <<
                  std::setw(8) << cp_iter << "|" <<
                  std::setw(8) << cp_flag << "|" <<
                  std::setw(10) << n_ineq << "|" <<
                  std::setw(12) << lb << "|" <<
                  std::setw(11) << glb << "*|" <<
                  std::setw(8) << set << "|" <<
                  std::setw(8) << i << " " <<
                  std::setw(8) << j << "|" <<
                  std::setw(13) << node_gap << "|" <<
                  std::setw(13) << gap << "|" <<
                  std::setw(8) << open << "|"
                  << std::endl;

    }
}