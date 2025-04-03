//
// Created by Rebecca Bekker on 3/26/25.
//

#ifndef CELL_H
#define CELL_H


#include <array>

class Cell {
public:
    Cell(std::array<double,2> loc, int type);
    std::array<double,2> x;
    int state;
};
#endif //CELL_H
