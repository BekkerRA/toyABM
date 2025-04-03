//
// Created by Rebecca Bekker on 3/26/25.
//
#include <Cell.h>
#include <_xlocale.h>

Cell::Cell(std::array<double, 2> loc, int cell_state) {
    x = loc;
    state = cell_state;
}

