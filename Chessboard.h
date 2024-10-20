// Andrei Shchapaniak
// 27.03.2024

#ifndef NI_PDP_CHESSBOARD_H
#define NI_PDP_CHESSBOARD_H

#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include <fstream>
#include <memory>

struct Point {
    int x, y;
};

struct Rectangle {
    Point p1, p2;
};

struct Knight {
    char color;
    int pos;
    bool active;

    Knight(char c, int p, bool a) : color(c), pos(p), active(a) {}

    bool operator<(const Knight& other) const {
        if (color != other.color) return color < other.color;
        if (pos != other.pos) return pos < other.pos;
        return active < other.active;
    }

    bool operator==(const Knight& other) const {
        return active == other.active && pos == other.pos && color == other.color;
    }
};

struct KnightsColor {
    std::vector<Knight> knights;
    std::vector<int> expected;
};

typedef struct Board {
    std::vector<char> board;
    std::map<char, KnightsColor> knightsc;
    char current_color;
    std::vector<char> colors;
    int active_knights;
    int n;
    int m;
} board_st;

class Chessboard {
public:
    std::vector<char> board;
    std::map<char, KnightsColor> knightsc;
    char current_color;
    std::vector<char> colors;

    explicit Chessboard(const std::string& filePath) : max_depth(0), current_color(1), colors({1 , 2}) {
        std::ifstream inputFile(filePath);
        if (!inputFile.is_open())
            throw std::runtime_error("Cannot open input file");

        inputFile >> m >> n >> k >> k;
        inputFile >> W.p1.x >> W.p1.y >> W.p2.x >> W.p2.y;
        inputFile >> B.p1.x >> B.p1.y >> B.p2.x >> B.p2.y;

        board.resize(n*m, 0);
        active_knights = 2*k;
        place_knights();
        calculate_max_depth();
    }

    explicit Chessboard(
        std::vector<char> board,
        std::map<char, KnightsColor> knightsc,
        char current_color,
        std::vector<char> colors,
        int active_knights,
        int n,
        int m
    ) : board(std::move(board)), knightsc(std::move(knightsc)),
        current_color(current_color), colors(std::move(colors)),
        active_knights(active_knights), n(n), m(m)
    {}

    explicit Chessboard() {}

    void set_active_knights(int value) {
        active_knights = value;
    }

    ~Chessboard() = default;

    [[nodiscard]] int get_m() const {
        return m;
    }

    [[nodiscard]] int get_n() const {
        return n;
    }

    [[nodiscard]] int get_active_knights() const {
        return active_knights;
    }

    [[nodiscard]] int get_max_depth() const {
        return max_depth;
    }

private:
    int m{}, n{}, k{}, max_depth{}, active_knights;
    Rectangle W{}, B{};

    void process_rectangle(Rectangle& rt, char color) {
        for (int x = rt.p1.x; x <= rt.p2.x; ++x) {
            for (int y = rt.p1.y; y <= rt.p2.y; ++y) {
                int pos = x * n + y;
                Knight knight = Knight(color, pos, true); // Use smart pointer

                board[pos] = color;
                knightsc[color].knights.push_back(knight);
                knightsc[3 - static_cast<int>(color)].expected.push_back(pos); // More explicit casting
            }
        }
    }

    void place_knights() {
        process_rectangle(W, 1);
        process_rectangle(B, 2);
    }

    [[nodiscard]] int calculate_distance(const Knight& knight, const std::vector<int>& final) const {
        int tmp_max_depth = 0;
        int startX, startY, finalX, finalY, distance, estimated;

        for (const int& fpos : final) {
            startX = knight.pos / n;
            startY = knight.pos % n;
            finalX = fpos / n;
            finalY = fpos % n;

            distance = std::abs(finalX - startX) + std::abs(finalY - startY);
            estimated = (distance + 2) / 3;

            tmp_max_depth = std::max(tmp_max_depth, estimated);
        }

        return tmp_max_depth + 1;
    }

    void calculate_max_depth() {
        for (const Knight& knight : knightsc.at(1).knights) {
            max_depth += calculate_distance(knight, knightsc.at(1).expected);
        }

        for (const Knight& knight : knightsc.at(2).knights) {
            max_depth += calculate_distance(knight, knightsc.at(2).expected);
        }
    }
};

#endif //NI_PDP_CHESSBOARD_H
