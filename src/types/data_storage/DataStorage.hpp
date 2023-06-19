/*
 * DataStorage.hpp
 *
 * This class is created to generic interface to store data.
 *  Created on: 28 Aug 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

#include <tuple>
#include <vector>

template <typename... Funcs>
class DataStorage
{
public:
  DataStorage(Funcs... funcs) : funcs_(std::make_tuple(funcs...)) {}

  template <typename T>
  void save(const T &su)
  {
    (std::apply([&](auto func) { data_.push_back(func(su)); }, funcs_), ...);
  }

private:
  std::tuple<Funcs...> funcs_;
  std::vector<double> data_;
};
