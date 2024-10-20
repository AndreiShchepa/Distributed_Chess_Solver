// Andrei Shchapaniak
// 27.03.2024

#ifndef NI_PDP_SOLVER_H
#define NI_PDP_SOLVER_H

#include <list>
#include <algorithm>
#include <utility>
#include <limits>
#include "Chessboard.h"
#include "State.h"
#include <mpi.h>
#include <omp.h>
#include <queue>
#include <thread>

#define TAG_DO_TASK 1
#define TAG_TASK_DONE 2
#define TAG_KILL 3
#define TAG_NEW_MAX_DEPTH 4

typedef struct PRet {
    int res_depth;
    std::list<std::vector<char>> solution;
} pret_t;

class Solver {
private:
    int max_depth, new_max_depth;
    bool solved;
    MPI_Request max_depth_request;
    std::vector<MPI_Request> requests;
    bool receiving = false;
    std::list<std::vector<char>> final_solution;

    [[nodiscard]] static int calc_dst(const int pos, const std::vector<int>& final, int n) {
        int tmp_min_depth = std::numeric_limits<int>::max();
        int startX = pos / n;
        int startY = pos % n;

        for (int fpos : final) {
            int finalX = fpos / n;
            int finalY = fpos % n;

            int distance = std::abs(finalX - startX) + std::abs(finalY - startY);
            int estimated = (distance + 2) / 3; // Assuming this estimates knight moves

            tmp_min_depth = std::min(tmp_min_depth, estimated);
        }

        return tmp_min_depth;
    }

    static int calculate_min_distance(const Chessboard& board) {
        int min_depth = 0;

        auto add_distances = [&min_depth, &board](const std::vector<Knight>& knights, const std::vector<int>& expected_positions) {
            for (const Knight& knight : knights)
                if (knight.active) min_depth += calc_dst(knight.pos, expected_positions, board.get_n());
        };

        if (board.knightsc.find(1) != board.knightsc.end())
            add_distances(board.knightsc.at(1).knights, board.knightsc.at(1).expected);

        if (board.knightsc.find(2) != board.knightsc.end())
            add_distances(board.knightsc.at(2).knights, board.knightsc.at(2).expected);

        return min_depth;
    }

    static bool is_in_final_state(Knight knight, const std::vector<int>& f) {
        return std::any_of(f.begin(), f.end(), [&knight](int fpos) {
            return knight.pos == fpos;
        });
    }

    static void generate_legal_moves(const Knight knight, std::vector<int>& moves, const Chessboard& board) {
        moves.clear();
        moves.reserve(8);

        const int x = knight.pos / board.get_n();
        const int y = knight.pos % board.get_n();

        constexpr std::pair<int, int> moveOffsets[] = {
        {2, 1}, {2, -1}, {-2, 1}, {-2, -1},
        {1, 2}, {1, -2}, {-1, 2}, {-1, -2}
        };

        for (const auto& [dx, dy] : moveOffsets) {
            const int newX = x + dx;
            const int newY = y + dy;

            if (newX >= 0 && newX < board.get_m() && newY >= 0 && newY < board.get_n()) {
                const int new_pos = newX * board.get_n() + newY;
                if (board.board[new_pos] == 0)
                    moves.push_back(new_pos);
            }
        }
    }

    static void update_colors(Chessboard* chessboard) {
        bool hasWhite = false, hasBlack = false;

        for (const auto& color_knights_pair : chessboard->knightsc) {
            for (const Knight& knight : color_knights_pair.second.knights) {
                if (!knight.active) continue;

                if (color_knights_pair.first == 1) hasWhite = true;
                else if (color_knights_pair.first == 2) hasBlack = true;

                if (hasWhite && hasBlack) break;
            }
        }

        if (hasWhite && hasBlack) return;

        chessboard->colors.clear();
        if (hasWhite) chessboard->colors.push_back(1);
        if (hasBlack) chessboard->colors.push_back(2);
    }

    static std::vector<std::pair<int, int>> sort_kmoves(Chessboard& chessboard, char color) {
        auto& knights = chessboard.knightsc.at(color).knights;
        auto& expectedPositions = chessboard.knightsc.at(color).expected;
        int boardSize = chessboard.get_n();

        std::vector<std::tuple<int, int, int>> distances; // <distance, index, move>
        for (int i = 0; i < knights.size(); i++) {
            if (!knights[i].active) continue;

            std::vector<int> moves;
            generate_legal_moves(knights[i], moves, chessboard);

            for (int move : moves)
                distances.emplace_back(calc_dst(move, expectedPositions, boardSize), i, move);
        }

        std::sort(distances.begin(), distances.end(), [](const auto& a, const auto& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        std::vector<std::pair<int, int>> sorted_pairs;
        sorted_pairs.reserve(distances.size());
        for (const auto& [distance, index, move] : distances)
            sorted_pairs.emplace_back(index, move);

        return sorted_pairs;
    }

    void check_final_solution(const int& depth, std::list<std::vector<char>>& solution_moves, int rank) {
        if (depth < max_depth) {
            #pragma omp critical
            if (depth < max_depth) {
                solved = true;
                max_depth = depth;
                final_solution = std::move(solution_moves);

                int proc_count;
                MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

                for (int i = 1; i < proc_count; i++) {
                    if (i != rank) {
                        MPI_Request new_request;
                        MPI_Isend(&depth, 1, MPI_INT, i, TAG_NEW_MAX_DEPTH, MPI_COMM_WORLD, &new_request);
                        requests.push_back(new_request);
                    }
                }
            }
        }
    }

    void start_receiving_for_max_depth() {
        if (!receiving) {
            MPI_Irecv(&new_max_depth, 1, MPI_INT, MPI_ANY_SOURCE, TAG_NEW_MAX_DEPTH, MPI_COMM_WORLD, &max_depth_request);
            receiving = true;
        }
    }

    void check_and_update_max_depth() {
        int flag = 0;
        MPI_Status status;
        MPI_Test(&max_depth_request, &flag, &status);
        if (flag) {
            MPI_Wait(&max_depth_request, MPI_STATUS_IGNORE);
            if (new_max_depth < max_depth)
                max_depth = new_max_depth;

            receiving = false;
        }
    }

    void dfs(State state, int rank) {
        start_receiving_for_max_depth();
        check_and_update_max_depth();

        if (state.depth >= 0.35 * max_depth) {
            no_parallel_dfs(State(state.depth, state.board, state.moves), rank);
            return;
        }

        for (std::pair<int, int> p : sort_kmoves(state.board, state.board.current_color)) {
            #pragma omp task firstprivate(state, p, rank) default(none)
            {
                Chessboard new_board = state.board;
                auto& knights = new_board.knightsc.at(new_board.current_color).knights;

                int prev_pos = knights[p.first].pos;
                new_board.board[prev_pos] = 0;
                new_board.board[p.second] = knights[p.first].color;
                knights[p.first].pos = p.second;

                update_colors(&new_board);

                if (is_in_final_state(knights[p.first], new_board.knightsc.at(new_board.current_color).expected)) {
                    knights[p.first].active = false;
                    new_board.set_active_knights(new_board.get_active_knights() - 1);
                }

                state.moves.push_back(new_board.board);

                if (new_board.colors.size() > 1)
                    new_board.current_color = state.board.current_color == 1 ? 2 : 1;

                dfs(State(state.depth + 1, new_board, state.moves), rank);
                state.moves.pop_back();
            }
        }
    }

    void no_parallel_dfs(State state, int rank) {
        start_receiving_for_max_depth();
        check_and_update_max_depth();

        if (state.board.get_active_knights() == 0) {
            check_final_solution(state.depth, state.moves, rank);
            return;
        }

        if (state.depth >= max_depth || state.depth + calculate_min_distance(state.board) + 1 >= max_depth)
            return;

        for (std::pair<int, int> p : sort_kmoves(state.board, state.board.current_color)) {
            Chessboard new_board = state.board;
            auto& knights = new_board.knightsc.at(new_board.current_color).knights;

            int prev_pos = knights[p.first].pos;
            new_board.board[prev_pos] = 0;
            new_board.board[p.second] = knights[p.first].color;
            knights[p.first].pos = p.second;

            update_colors(&new_board);

            if (is_in_final_state(knights[p.first], new_board.knightsc.at(new_board.current_color).expected)) {
                knights[p.first].active = false;
                new_board.set_active_knights(new_board.get_active_knights() - 1);
            }

            state.moves.push_back(new_board.board);

            if (new_board.colors.size() > 1)
                new_board.current_color = state.board.current_color == 1 ? 2 : 1;

            no_parallel_dfs(State(state.depth + 1, new_board, state.moves), rank);

            state.moves.pop_back();
        }
    }

public:
    explicit Solver(const Chessboard& board)
            : max_depth(board.get_max_depth()), solved(false) {}

    void print_solution(int n) {
        char num_let[3] = {'.', 'W', 'B'};
        if (!solved) {
            printf("Solution was not found!\n");
            return;
        }

        printf("%d\n\n", max_depth);

        for (std::vector<char>& state: final_solution) {
            for (int i = 0; i < state.size(); i++) {
                printf("%c ", num_let[state[i]]);
                if ((i+1) % n == 0) printf("\n");
            }
            printf("\n");
        }
    }

    [[nodiscard]] std::queue<State> generate_states(int max_count, const Chessboard& board, const std::list<std::vector<char>>& moves) const {
        std::vector<State> states_v;
        std::queue<State> states_q;
        int depth = 0;
        states_q.emplace(depth, board, moves);

        while (states_q.size() < max_count) {
            State curr_s = states_q.front();
            if (board.get_active_knights() == 0)
                break;

            if (depth >= max_depth/2)
                break;

            states_q.pop();
            for (std::pair<int, int> p : sort_kmoves(curr_s.board, curr_s.board.current_color)) {
                Chessboard new_board = curr_s.board;
                auto& knights = new_board.knightsc.at(new_board.current_color).knights;

                int prev_pos = knights[p.first].pos;
                new_board.board[prev_pos] = 0;
                new_board.board[p.second] = knights[p.first].color;
                knights[p.first].pos = p.second;

                update_colors(&new_board);

                if (is_in_final_state(knights[p.first], new_board.knightsc.at(new_board.current_color).expected)) {
                    knights[p.first].active = false;
                    new_board.set_active_knights(new_board.get_active_knights() - 1);
                }

                curr_s.moves.push_back(new_board.board);

                if (new_board.colors.size() > 1)
                    new_board.current_color = curr_s.board.current_color == 1 ? 2 : 1;

                states_q.emplace( curr_s.depth + 1, new_board, curr_s.moves);
                curr_s.moves.pop_back();
            }
        }

        return states_q;
    }

    static std::vector<int> serialize(const pret_t& data) {
        std::vector<int> serialized;

        serialized.push_back(data.res_depth);
        serialized.push_back(data.solution.size());
        for (const auto& move : data.solution) {
            serialized.push_back(move.size());
            for (char c : move) {
                serialized.push_back(static_cast<int>(c));
            }
        }

        return serialized;
    }

    static pret_t deserialize(const std::vector<int>& data) {
        pret_t ret_obj;
        size_t idx = 0;
        ret_obj.res_depth = data[idx++];

        size_t movesSize = data[idx++];
        ret_obj.solution.clear();
        for (size_t i = 0; i < movesSize; ++i) {
            size_t moveSize = data[idx++];
            std::vector<char> move;
            move.reserve(moveSize);
            for (size_t j = 0; j < moveSize; ++j) {
                move.push_back(static_cast<char>(data[idx++]));
            }
            ret_obj.solution.push_back(move);
        }

        return ret_obj;
    }

    void master(const int& proc_count, const Chessboard& chessboard, int task_size, int res_size) {
        int worker_count = proc_count - 1;
        std::list<std::vector<char>> solution_moves = {};

        solution_moves.emplace_back(chessboard.board);
        std::queue<State> states = generate_states(proc_count * 1, chessboard, solution_moves);

        for (int i = 1; i <= worker_count && !states.empty(); i++) {
            State task = states.front();
            std::vector<int> data_buffer = task.serialize();
            data_buffer.resize(task_size);
            states.pop();

            MPI_Send(data_buffer.data(), task_size, MPI_INT, i, TAG_DO_TASK, MPI_COMM_WORLD);
        }

        while (worker_count > 0) {
            pret_t result;
            std::vector<int> byte_result(res_size);
            MPI_Status status;
            MPI_Recv(byte_result.data(), res_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            result = deserialize(byte_result);

            if (status.MPI_TAG == TAG_TASK_DONE) {
                if (result.solution.size() > 0 && (result.res_depth < max_depth || final_solution.size() > result.solution.size())) {
                    max_depth = result.res_depth;
                    final_solution = std::move(result.solution);
                    solved = true;
                }

                if (states.empty()) {
                    int nothing = 0;
                    MPI_Send(&nothing, 1, MPI_INT, status.MPI_SOURCE, TAG_KILL, MPI_COMM_WORLD);
                    worker_count--;
                    continue;
                }

                State task = states.front();
                states.pop();
                std::vector<int> data_buffer = task.serialize();
                data_buffer.resize(task_size);
                MPI_Send(data_buffer.data(), task_size, MPI_INT, status.MPI_SOURCE, TAG_DO_TASK,MPI_COMM_WORLD);
            }
        }

        print_solution(chessboard.get_n());
    }

    void slave(int rank, int task_size, int res_size) {
        while (true) {
            std::vector<int> byte_result(task_size);
            MPI_Status status;
            MPI_Recv(byte_result.data(), task_size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == TAG_KILL) {
                if (receiving) {
                    MPI_Cancel(&max_depth_request);
                    MPI_Request_free(&max_depth_request);

                    for (auto& req : requests) {
                        MPI_Cancel(&req);
                        MPI_Request_free(&req);
                    }
                }

                break;
            }

            if (status.MPI_TAG == TAG_DO_TASK) {
                pdata_t result = State().deserialize(byte_result);

                dfs(State(
                        result.depth,
                  Chessboard(
                          result.board.board, result.board.knightsc,
                          result.board.current_color, result.board.colors,
                          result.board.active_knights,
                          result.board.n, result.board.m),
           result.moves), rank);

                std::vector<int> send_res = serialize({this->max_depth, this->final_solution});
                send_res.resize(res_size);
                MPI_Send(send_res.data(), res_size, MPI_INT, 0, TAG_TASK_DONE,MPI_COMM_WORLD);
            }
        }
    }

    void start(const Chessboard& chessboard) {
        int rank, proc_count;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

        int task_size = 12 + 3 * chessboard.get_n() * chessboard.get_m() + 6 * chessboard.get_active_knights();
        int res_size = 2 + chessboard.get_n() * chessboard.get_m() * max_depth;

        if (rank == 0) {
            master(proc_count, chessboard, task_size, res_size);
        } else {
            slave(rank, task_size, res_size);
        }
    }
};

#endif //NI_PDP_SOLVER_H
