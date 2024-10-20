// Andrei Shchapaniak
// 27.03.2024

#ifndef NI_PDP_STATE_H
#define NI_PDP_STATE_H

#include <utility>
#include <vector>
#include <list>
#include "Chessboard.h"

typedef struct PData {
    int depth;
    board_st board;
    std::list<std::vector<char>> moves;
} pdata_t;


class State {
public:
    int depth{};
    Chessboard board;
    std::list<std::vector<char>> moves;

    explicit State(int depth, const Chessboard& board, std::list<std::vector<char>> solution_moves) :
        depth(depth),
        board(board),
        moves(std::move(solution_moves)) {}

    explicit State() = default;

    std::vector<int> serialize() {
        std::vector<int> serialized;

        serialized.push_back(board.board.size());
        for (char c : board.board) {
            serialized.push_back(static_cast<int>(c));
        }

        serialized.push_back(board.knightsc.size());
        for (const auto& pair : board.knightsc) {
            serialized.push_back(static_cast<int>(pair.first)); // Key (char to int)
            const auto& knightsColor = pair.second;

            serialized.push_back(knightsColor.knights.size());
            for (const auto& knight : knightsColor.knights) {
                serialized.push_back(static_cast<int>(knight.color));
                serialized.push_back(knight.pos);
                serialized.push_back(knight.active ? 1 : 0);
            }

            serialized.push_back(knightsColor.expected.size());
            for (int expected : knightsColor.expected) {
                serialized.push_back(expected);
            }
        }

        serialized.push_back(static_cast<int>(board.current_color));
        serialized.push_back(board.colors.size());
        for (char c : board.colors) {
            serialized.push_back(static_cast<int>(c));
        }

        serialized.push_back(board.get_active_knights());
        serialized.push_back(board.get_n());
        serialized.push_back(board.get_m());

        serialized.push_back(depth);
        serialized.push_back(moves.size());
        for (const auto& move : moves) {
            serialized.push_back(move.size());
            for (char c : move) {
                serialized.push_back(static_cast<int>(c));
            }
        }

        return serialized;
    }

    pdata_t deserialize(const std::vector<int>& serialized) {
        pdata_t ret_data;
        size_t idx = 0;

        size_t boardSize = serialized[idx++];
        ret_data.board.board.clear();
        ret_data.board.board.reserve(boardSize);
        for (size_t i = 0; i < boardSize; ++i) {
            ret_data.board.board.push_back(static_cast<char>(serialized[idx++]));
        }

        size_t knightscSize = serialized[idx++];
        ret_data.board.knightsc.clear();
        for (size_t i = 0; i < knightscSize; ++i) {
            char key = static_cast<char>(serialized[idx++]);
            KnightsColor knightsColor;

            size_t knightsSize = serialized[idx++];
            for (size_t j = 0; j < knightsSize; ++j) {
                char color = static_cast<char>(serialized[idx++]);
                int pos = serialized[idx++];
                bool active = serialized[idx++] == 1;
                knightsColor.knights.emplace_back(color, pos, active);
            }

            size_t expectedSize = serialized[idx++];
            for (size_t j = 0; j < expectedSize; ++j) {
                knightsColor.expected.push_back(serialized[idx++]);
            }

            ret_data.board.knightsc[key] = knightsColor;
        }

        ret_data.board.current_color = static_cast<char>(serialized[idx++]);

        size_t colorsSize = serialized[idx++];
        ret_data.board.colors.clear();
        ret_data.board.colors.reserve(colorsSize);
        for (size_t i = 0; i < colorsSize; ++i) {
            ret_data.board.colors.push_back(static_cast<char>(serialized[idx++]));
        }

        ret_data.board.active_knights = serialized[idx++];
        ret_data.board.n = serialized[idx++];
        ret_data.board.m = serialized[idx++];

        ret_data.depth = serialized[idx++];

        size_t movesSize = serialized[idx++];
        ret_data.moves.clear();
        for (size_t i = 0; i < movesSize; ++i) {
            size_t moveSize = serialized[idx++];
            std::vector<char> move;
            move.reserve(moveSize);
            for (size_t j = 0; j < moveSize; ++j) {
                move.push_back(static_cast<char>(serialized[idx++]));
            }
            ret_data.moves.push_back(move);
        }

        return ret_data;
    }
};

#endif //NI_PDP_STATE_H
