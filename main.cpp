// Andrei Shchapaniak
// 27.03.2024

#include <chrono>
#include <iostream>
#include "Chessboard.h"
#include "Solver.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>\n";
        return 1;
    }
       
    omp_set_num_threads(32);
    int procNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &procNum);

    Chessboard board = Chessboard(argv[1]);
    Solver solver = Solver(board);

    if (procNum != 0) {
        solver.start(board);
        MPI_Finalize();
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();
    solver.start(board);
    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = (finish - start);
    std::cout << "Duration: ";
    if (duration.count() <= 1000.0) {
        auto durationInMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        std::cout << durationInMilliSeconds << " ms" << std::endl;
    } else {
        auto durationInSeconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
        std::cout << durationInSeconds << " s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

